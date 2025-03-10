/** nelmin.go
Nelder-Mead simplex minimization of a nonlinear (multivariate) function.

This code has been adapted from the C-coded nelmin.c which was
adapted from the Fortran-coded nelmin.f which was, in turn, adapted
from the papers:

    J.A. Nelder and R. Mead (1965)
    A simplex method for function minimization.
    Computer Journal, Volume 7, pp 308-313.

    R. O'Neill (1971)
    Algorithm AS47. Function minimization using a simplex algorithm.
    Applied Statistics, Volume 20, pp 338-345.

For the 2020 update, we make some of the stepping closer to
the description given in the paper:

    Donghoon Lee and Matthew Wiswall (2007)
    A parallel implementation of the simplec function minimization routine.
    Computational Economics 30:171-187

Some examples are in:

   D.M. Olsson and L.S. Nelson (1975)
   The Nelder-Mead Simplex procedure for function minimization.
   Technometrics, Volume 17 No. 1, pp 45-51.

For a fairly recent and popular incarnation of this minimizer,
see the amoeba function in the famous "Numerical Recipes" text.
The programming interface is via the minimize() function; see below.

Author:
   PA Jacobs, School of Engineering, The University of Queensland

Version:
   2004-01-07 Python2 flavour for the cfcfd project.
   2020-06-22 Port to Python3 for the DGD project.
   2020-06-25 Concurrent evaluation of the candidate points.
   2021-06-07 Dan Smith added option to read the initial simplex.
   2024-01-15 Golang version
*/

package nelmin

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"sort"
	"github.com/pajacobs-ghub/nm/array"
)

//-----------------------------------------------------------------------------

type Vertex struct {
	X *array.Vector
	F float64
}

func NewVertex(n int) Vertex {
	return Vertex{X: array.NewVector(n), F: 0.0}
}

func (v Vertex) String() string {
	// We write a JSON compatible representation
	// so that we can use it when writing the full simplex.
	// It is readable enough for standard printing.
	return fmt.Sprintf("{%q:%s, %q:%g}", "x", v.X.String(), "f", v.F)
}

func approxEquals(a float64, b float64, tol float64) bool {
	// Relative comparison for large numbers, absolute comparison for small numbers
	return math.Abs(a-b)/(0.5*(math.Abs(a)+math.Abs(b)+1.0)) <= tol
}

func (v Vertex) ApproxEquals(other Vertex, tol float64) bool {
	return v.X.ApproxEquals(other.X, tol) && approxEquals(v.F, other.F, tol)
}

// We want to compute a centroid of a subset of vertices in the simplex,
// so leave the last p vertices out of the calculation.
func Centroid(va []Vertex, p int) (*Vertex, error) {
	nv := len(va)
	if nv == 0 {
		return nil, errors.New("No vertices in array.")
	}
	// Number of points with which to compute the centroid.
	k := nv - p
	if k <= 0 {
		return nil, errors.New("Not enough vertices remaining.")
	}
	n := len(va[0].X.Data)
	c := NewVertex(n)
	for i := 0; i < k; i++ {
		c.X.Add(c.X, va[i].X)
		c.F += va[i].F
	}
	s := 1.0/float64(k)
	c.X.Scale(s)
	c.F *= s
	return &c, nil
}

//-----------------------------------------------------------------------------

func MakeSimplexAboutPoint(
	f func([]float64) float64,
	x0 []float64,
	dx []float64) ([]Vertex, int, error) {
	n := len(x0)
	nfe := 0
	if n == 0 {
		return nil, nfe, errors.New("Zero number of parameters.")
	}
	if n != len(dx) {
		return nil, nfe, errors.New("len(dx) did not match len(x)")
	}
	anyZero := false
	for i := 0; i < n; i++ {
		if dx[i] == 0.0 {
			anyZero = true
		}
	}
	if anyZero {
		return nil, nfe, errors.New("One or more zero value in dx.")
	}
	// TODO Option to do the function evaluations in parallel.
	fx0 := f(x0)
	nfe += 1
	smplx := []Vertex{Vertex{array.NewVectorFromArray(x0), fx0}}
	for i := 0; i < n; i++ {
		x1 := make([]float64, n)
		for j := 0; j < n; j++ {
			x1[j] = x0[j]
		}
		x1[i] += dx[i]
		fx1 := f(x1)
		nfe += 1
		smplx = append(smplx, Vertex{array.NewVectorFromArray(x1), fx1})
	}
	sortSimplex(smplx)
	return smplx, nfe, nil
}

func sortSimplex(smplx []Vertex) {
	if len(smplx) < 2 {
		return
	}
	sort.Slice(smplx, func(i int, j int) bool {
		return smplx[i].F < smplx[j].F
	})
	return
}

func fStats(smplx []Vertex) (float64, float64, error) {
	// Returns mean and standard-deviation for the function values.
	nv := len(smplx)
	if nv == 0 {
		return 0.0, 0.0, errors.New("No points in simplex")
	}
	mean := 0.0
	for i := 0; i < nv; i++ {
		mean += smplx[i].F
	}
	mean /= float64(nv)
	variance := 0.0
	for i := 0; i < nv; i++ {
		dF := smplx[i].F - mean
		variance += dF * dF
	}
	variance /= float64(nv)
	return mean, math.Sqrt(variance), nil
}

func SimplexToJSON(smplx []Vertex) string {
	var b bytes.Buffer
	b.WriteString(fmt.Sprintf("{%q:%d, %q:%s}", "n", len(smplx)-1, "vertices", smplx))
	return b.String()
}

func SimplexFromJSON(str string) ([]Vertex, error) {
	smplx := []Vertex{}
	// TODO
	return smplx, nil
}

//-----------------------------------------------------------------------------

type Minimizer struct {
	F                func(x []float64) float64 // Client-supplied objective function.
	Vertices         []Vertex                  // The simplex is N+1 Vertices, where N is len(x).
	P                int                       // Number of points to be replaced in parallel.
	Steps            int                       // Steps between convergence checks.
	NFEvaluationsMax int                       // Limit function evaluations.
	NFEvaluations    int
	Nrestarts        int
	Kreflect         float64
	Kextend          float64
	Kcontract        float64
	Tol              float64
}

func NewMinimizer(f func([]float64) float64) *Minimizer {
	m := Minimizer{F: f,
		Vertices:         nil,
		P:                1,
		Steps:            20,
		NFEvaluationsMax: 300,
		NFEvaluations:    0,
		Nrestarts:        0,
		Kreflect:         1.0,
		Kextend:          2.0,
		Kcontract:        0.5,
		Tol:              1.0e-6}
	return &m
}

func (m *Minimizer) String() string {
	// Returns a JSON compatible string.
	return fmt.Sprintf("{%q:%p, %q:%s, %q:%d, %q:%d, %q:%d, %q:%d, %q:%d, %q:%g, %q:%g, %q:%g, %q:%g}",
		"fun", m.F, "vertices", m.Vertices, "p", m.P, "steps", m.Steps,
		"nfemax", m.NFEvaluationsMax, "nfe", m.NFEvaluations, "nrestarts", m.Nrestarts,
		"reflect", m.Kreflect, "extend", m.Kextend, "contract", m.Kcontract,
		"tol", m.Tol)
}

func (m *Minimizer) replaceVertex(i int, xMid *array.Vector) (bool, int) {
	// Try to replace the specified i vertex with a better point,
	// returning a flag to indicate if successful.
	//
	// Note that we may want to call this method concurrently
	// so that we can replace m.P points in parallel, being careful
	// that we don't try to replace the same point in more than one thread.
	// Also, note that the objective function calls will need to be truly
	// independent to make this work reliably.
	nfe := 0
	// Assuming a sorted array, 0 is the best point (minimum value of F).
	fMin := m.Vertices[0].F
	// Assume that Vertex[i] is a high value of the objective fn.
	xHigh := m.Vertices[i].X
	fHigh := m.Vertices[i].F
	// First, try moving away from worst point by reflection through centroid.
	n := len(xHigh.Data)
	xRefl := array.NewVector(n)
	xRefl.Blend(xMid, xHigh, (1.0+m.Kreflect), -m.Kreflect)
	fRefl := m.F(xRefl.Data)
	nfe += 1
	if fRefl < fMin {
		// The reflection through the centroid is good,
		// try to extend in the same direction.
		xExt := array.NewVector(n)
		xExt.Blend(xMid, xRefl, (1.0-m.Kextend), m.Kextend)
		fExt := m.F(xExt.Data)
		nfe += 1
		if fExt < fRefl {
			// Keep the extension because it's best.
			m.Vertices[i] = Vertex{xExt, fExt}
			return true, nfe
		} else {
			// Settle for the original reflection.
			m.Vertices[i] = Vertex{xRefl, fRefl}
			return true, nfe
		}
	} else {
		// The reflection is not going in the right direction, it seems.
		// See how many vertices are worse than the reflected point.
		count := 0
		for j := 0; j < n+1; j++ {
			if m.Vertices[j].F > fRefl {
				count += 1
			}
		}
		if count <= 1 {
			// Not too many points are higher than the original reflection.
			// Try a contraction on the reflection-side of the centroid.
			xCon := array.NewVector(n)
			xCon.Blend(xMid, xHigh, (1.0-m.Kcontract), m.Kcontract)
			fCon := m.F(xCon.Data)
			nfe += 1
			if fCon < fHigh {
				// At least we haven't gone uphill; accept.
				m.Vertices[i] = Vertex{xCon, fCon}
				return true, nfe
			}
		} else {
			// Retain the original reflection because there are many
			// original vertices with higher values of the objective function
			// and it will be good to have some change to the simplex.
			m.Vertices[i] = Vertex{xRefl, fRefl}
			return true, nfe
		}
	}
	// If we arrive here, we have not replaced the highest point.
	return false, nfe
} // end replaceVertex()

func (m *Minimizer) contractAboutBestPoint() {
	// Assuming a sorted array, 0 is the best point (minimum value of F).
	xMin := m.Vertices[0].X
	// Move all other simplex vertices to half-way between their current point
	// and the best point.
	// TODO Option to do the function evaluations in parallel.
	nv := len(m.Vertices)
	for i := 1; i < nv; i++ {
		m.Vertices[i].X.Blend(xMin, m.Vertices[i].X, 0.5, 0.5)
		m.Vertices[i].F = m.F(m.Vertices[i].X.Data)
	}
	m.NFEvaluations += nv-1
	return
}

func (m *Minimizer) TakeSteps(nsteps int) error {
	// Take some steps, updating the simplex.
	// On return, the best point is m.Vertices[0].
	nv := len(m.Vertices)
	for step := 0; step < nsteps; step++ {
		// Compute the centroid of the points that we are not replacing.
		vMid, err := Centroid(m.Vertices, m.P)
		if err != nil {
			return fmt.Errorf("Error while computing centroid: %s", err)
		}
		// Try to replace the P worst points by generating new points
		// about the current centroid (vMid).
		// TODO Option to do the function evaluations in parallel.
		anySuccess := false
		for i := 0; i < m.P; i++ {
			success, nfe := m.replaceVertex(nv-1-i, vMid.X)
			if success {
				anySuccess = true
			}
			m.NFEvaluations += nfe
		}
		if !anySuccess {
			// Did not improve any of the worst points.
			m.contractAboutBestPoint()
			m.Nrestarts += 1
		}
		sortSimplex(m.Vertices)
	}
	return nil
}

func (m *Minimizer) MinimizeFromPoint(x []float64, dx []float64) error {
	var err error
	var nfe int
	m.Vertices, nfe, err = MakeSimplexAboutPoint(m.F, x, dx)
	if err != nil {
		return fmt.Errorf("Error while making initial simplex: %s", err)
	}
	m.NFEvaluations += nfe
	for m.NFEvaluations < m.NFEvaluationsMax {
		m.TakeSteps(m.Steps)
		_, sdev, err := fStats(m.Vertices)
		if err != nil {
			return fmt.Errorf("Error while computing function stats: %s", err)
		}
		if sdev < m.Tol {
			// Points within the simplex have similar function values,
			// and we deem this to be good enough to stop stepping.
			break
		}
	}
	return nil
}
