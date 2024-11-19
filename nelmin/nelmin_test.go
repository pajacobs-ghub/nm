/** nelmin_test.go
 * Try out the Nelder-Mead simplex optimizer.
 *
 * Peter J, 2024-Jan-15
 */

package nelmin

import (
	"fmt"
	"math"
	"testing"
)

func TestVertex(t *testing.T) {
	fmt.Println("Test Vertex Functions")
	v1 := Vertex{X: []float64{1.0, 2.0, 3.0}, F: 99.0}
	v2 := Vertex{[]float64{1.0, 2.0, 3.0}, 99.0}
	if !v1.ApproxEquals(v2, 1.0e-6) {
		t.Errorf("Initialization: Should be the same v1=%v, v2=%v", v1, v2)
	}
	v3 := NewVertex(3)
	if v1.ApproxEquals(v3, 1.0e-6) {
		t.Errorf("NewVertex: Should be different v1=%v, v3=%v", v1, v3)
	}
	v4, err := Centroid([]Vertex{v1, v2, v3}, 1)
	// v4, err := Centroid([]Vertex{}) // To exercise error path.
	if err != nil {
		t.Errorf("Centroid calc, got error: %s", err)
	} else {
		v5 := Vertex{X: []float64{1.0, 2.0, 3.0}, F: 99.0}
		if !v4.ApproxEquals(v5, 1.0e-6) {
			t.Errorf("Centroid calc, Should be the same v4=%v, v5=%v", v1, v2)
		}
	}
	v6 := NewVertex(3)
	v6, err = v6.Add(v1, 1.0)
	if !v6.ApproxEquals(v2, 1.0e-6) {
		t.Errorf("Add with scale: Should be the same v6=%v, v2=%v", v6, v2)
	}
}

func obj1(x []float64) float64 {
	n := len(x)
	s := 0.0
	for i := 0; i < n; i++ {
		s += (x[i] - 1.0) * (x[i] - 1.0)
	}
	return s
}

func TestSimplex(t *testing.T) {
	fmt.Println("Test Simplex functions")
	x0 := []float64{1.0, 2.0, 3.0}
	dx := []float64{0.1, 0.2, 0.3}
	smplx, nfe, err := MakeSimplexAboutPoint(obj1, x0, dx)
	if err != nil {
		t.Errorf("Failed to make simplex, err: %s", err)
	}
	// fmt.Println("smplx=", SimplexToJSON(smplx))
	// fmt.Println("nfe=", nfe)
	if nfe != 4 {
		t.Errorf("Simplex setup: nfe=%v should be 4", nfe)
	}
	var vmid *Vertex
	vmid, err = Centroid(smplx, 1)
	// fmt.Println("vmid=", vmid)
	vRef := Vertex{X: []float64{1.0333, 2.0667, 3.0}, F: 5.150}
	if !vmid.ApproxEquals(vRef, 1.0e-3) {
		t.Errorf("Simplex centroid: Should be same vmid=%v, vRef=%v", vmid, vRef)
	}
}

func TestMinimizer1(t *testing.T) {
	fmt.Println("Test Minimizer functions for simple quadratic objective")
	m := NewMinimizer(obj1)
	x := []float64{0.0, 0.0, 0.0}
	dx := []float64{0.1, 0.1, 0.1}
	err := m.MinimizeFromPoint(x, dx)
	if err != nil {
		t.Errorf("Failed to mimimize from point, err: %s", err)
	}
	// fmt.Printf("m=%s\n", m.String())
	if m.NFEvaluations != 106 {
		t.Errorf("Wrong number of function evaluations: m.nfe=%v should be 106", m.NFEvaluations)
	}
	if m.Nrestarts != 0 {
		t.Errorf("Wrong number of restarts: nrestarts=%v should be 0", m.Nrestarts)
	}
	vRef := Vertex{X: []float64{1.0, 1.0, 1.0}, F: 0.0}
	vMin := m.Vertices[0]
	if !vMin.ApproxEquals(vRef, 1.0e-3) {
		t.Errorf("Wrong minimum point: Should be same vMin=%v, vExpected=%v", vMin, vRef)
	}
}

func obj2(x []float64) float64 {
	// Example 3.3 from Olsson and Nelson
	x1, x2 := x[0], x[1] // Rename to match paper.
	if (x1*x1 + x2*x2) > 1.0 {
		return 1.0e38
	} else {
		yp := 53.69 + 7.26*x1 - 10.33*x2 + 7.22*x1*x1 + 6.43*x2*x2 + 11.36*x1*x2
		ys := 82.17 - 1.01*x1 - 8.61*x2 + 1.40*x1*x1 - 8.76*x2*x2 - 7.20*x1*x2
		return -yp + math.Abs(ys-87.8)
	}
}

func TestMinimizer2(t *testing.T) {
	fmt.Println("Example 3.3 from Olsson and Nelson")
	m := NewMinimizer(obj2)
	m.Tol = 1.0e-4
	x := []float64{0.0, 0.0}
	dx := []float64{0.5, 0.5}
	err := m.MinimizeFromPoint(x, dx)
	if err != nil {
		t.Errorf("Failed to mimimize from point, err: %s", err)
	}
	if m.NFEvaluations != 82 {
		t.Errorf("Wrong number of function evaluations: nfe=%v should be 82", m.NFEvaluations)
	}
	if m.Nrestarts != 0 {
		t.Errorf("Wrong number of restarts: nrestarts=%v should be 0", m.Nrestarts)
	}
	vRef := Vertex{X: []float64{0.811, -0.585}, F: -67.1}
	vMin := m.Vertices[0]
	if !vMin.ApproxEquals(vRef, 1.0e-3) {
		t.Errorf("Wrong minimum point, Should be same vMin=%v, vRef=%v", vMin, vRef)
	}
	// fmt.Printf("m=%s\n", m.String())
}

func obj3(z []float64) float64 {
	// Example 3.5 from Olsson and Nelson, least squares
	a1, a2, alpha1, alpha2 := z[0], z[1], z[2], z[3]
	x := []float64{0.25, 0.50, 1.00, 1.70, 2.00, 4.00}
	y := []float64{0.25, 0.40, 0.60, 0.58, 0.54, 0.27}
	s := 0.0
	for i := 0; i < len(x); i++ {
		t := x[i]
		eta := a1*math.Exp(alpha1*t) + a2*math.Exp(alpha2*t)
		r := y[i] - eta
		s += r * r
	}
	return s
}

func TestMinimizer3(t *testing.T) {
	fmt.Println("Example 3.5 from Olsson and Nelson, least squares")
	m := NewMinimizer(obj3)
	m.NFEvaluationsMax = 800
	m.Tol = 1.0e-9
	m.P = 2
	x := []float64{1.0, 1.0, -0.5, -2.5}
	dx := []float64{0.1, 0.1, 0.1, 0.1}
	err := m.MinimizeFromPoint(x, dx)
	if err != nil {
		t.Errorf("Failed to mimimize from point, err: %s", err)
	}
	if m.NFEvaluations != 495 {
		t.Errorf("Wrong number of function evaluations: nfe=%v should be 495", m.NFEvaluations)
	}
	if m.Nrestarts != 0 {
		t.Errorf("Wrong number of restarts: nrestarts=%v should be 0", m.Nrestarts)
	}
	vRef := Vertex{X: []float64{1.801, -1.842, -0.463, -1.205}, F: 0.0009}
	vMin := m.Vertices[0]
	if !vMin.ApproxEquals(vRef, 1.0e-3) {
		t.Errorf("Example 3.5, Should be the same vMin=%v, vRef=%v", vMin, vRef)
	}
	// fmt.Printf("m=%s\n", m.String())
}
