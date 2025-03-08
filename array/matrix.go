// matrix.go
// A package for working with augmented matrices for linear algebra.
// 
// The underlying matrix storage is an array of arrays, row-major convention,
// and we use the text-book description of elimination to find the inverse.
// Intended use is for small-scale exercises. For anything big, prefer gonum.
//
// Peter J. 2025-03-07
//

package array

import (
	"fmt"
	"errors"
	"math"
	"bytes"
)

type Matrix struct {
	Data [][]float64
}

func (v Matrix) IsEmpty() bool {
	return len(v.Data) == 0
}

func NewMatrix(nrows, ncols int) (*Matrix, error) {
	if nrows <= 0 {
		msg := fmt.Sprintf("Invalid value for nrows=%v", nrows)
		return nil, errors.New(msg)
	}
	z := Matrix{Data: make([][]float64, nrows)}
	if ncols <= 0 {
		msg := fmt.Sprintf("Invalid value for ncols=%v", ncols)
		return &z, errors.New(msg)
	}
	for i := 0; i < nrows; i++ {
		z.Data[i] = make([]float64, ncols)
	}
	return &z, nil
}

func NewMatrixFromArray(data [][]float64) (*Matrix, error) {
	z := Matrix{}
	nrows := len(data)
	if nrows == 0 {
		return &z, errors.New("Zero rows")
	}
	z.Data = make([][]float64, nrows)
	ncols0 := len(data[0])
	if ncols0 == 0 {
		return &z, errors.New("Zero columns")
	}
	for i := 0; i < nrows; i++ {
		ncols := len(data[i])
		if ncols != ncols0 {
			msg := fmt.Sprintf("Ragged rows: ncols0=%d ncols[%d]=%d", ncols0, i, ncols)
			return &z, errors.New(msg)
		}
		z.Data[i] = make([]float64, ncols0)
		for j := 0; j < ncols; j++ {
			z.Data[i][j] = data[i][j]
		}
	}
	return &z, nil
}

func (a *Matrix) String() string {
	var b bytes.Buffer
	nrows := len(a.Data)
	b.WriteString("[")
	for i := 0; i < nrows; i++ {
		ncols := len(a.Data[i])
		b.WriteString("[")
		for j := 0; j < ncols; j++ {
			b.WriteString(fmt.Sprintf("%g", a.Data[i][j]))
			if j+1 < ncols {
				b.WriteString(", ")
			}
		}
		b.WriteString("]")
		if i+1 < nrows {
			b.WriteString(", ")
		}
	}
	b.WriteString("]")
	return b.String()
}

func (a *Matrix) ApproxEquals(other *Matrix, tol float64) bool {
	nrows := len(a.Data)
	if nrows != len(other.Data) {
		return false
	}
	for i := 0; i < nrows; i++ {
		ncols := len(a.Data[i])
		if ncols != len(other.Data[i]) {
			return false
		}
		for j := 0; j < ncols; j++ {
			aa := a.Data[i][j]
			bb := other.Data[i][j]
	    	// Relative comparison for large numbers, absolute comparison for small numbers
			if math.Abs(aa-bb)/(0.5*(math.Abs(aa)+math.Abs(bb)+1.0)) > tol {
				return false
			}
		}
	}
	return true
}

func (a *Matrix) NormInf() float64 {
	norm := 0.0
	nrows := len(a.Data)
	if nrows == 0 {
		return norm
	}
	ncols := len(a.Data[0])
	if ncols == 0 {
		return norm
	}
	for i := 0; i < nrows; i++ {
		rowsum := 0.0
		for j := 0; j < ncols; j++ {
			rowsum += math.Abs(a.Data[i][j])
		}
		norm = math.Max(rowsum, norm)
	}
	return norm
}

var verySmallValue float64 = 1.0e-16

func SetVerySmallValue(v float64) {
	verySmallValue = v
	return
}

// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
// When computing an inverse, the incoming data is assumed to be c=[A|I].
func (c *Matrix) GaussJordanElimination() (*Matrix, error) {
	nrows := len(c.Data)
	if nrows == 0 {
		return c, errors.New("Empty Matrix")
	}
	ncols := len(c.Data[0])
	if ncols == 0 {
		return c, errors.New("Empty rows in Matrix")
	}
	if nrows == 1 {
		det := c.Data[0][0]
		if math.Abs(det) < verySmallValue {
			return c, errors.New(fmt.Sprintf("Singular with det=%v", det))
		}
		for j := 0; j < ncols; j++ {
			c.Data[0][j] /= det
		}
		return c, nil
	}
	for j := 0; j < nrows; j++ {
		// Select pivot, the largest magnitude in column j.
		p := j
		for i := j+1; i < nrows; i++ {
			if math.Abs(c.Data[i][j]) > math.Abs(c.Data[p][j]) {
				p = i
			}
		}
		if math.Abs(c.Data[p][j]) < verySmallValue {
			return c, errors.New(fmt.Sprintf("Singular with pivot=%v", c.Data[p][j]))
		}
		if p != j {
			// Swap rows, to get pivot onto the diagonal.
			c.Data[p], c.Data[j] = c.Data[j], c.Data[p]
		}
		// Scale row j to get unity on the diagonal.
		cjj := c.Data[j][j]
		for col := 0; col < ncols; col++ {
			c.Data[j][col] /= cjj
		}
		// Do the elimination to get zeros in all off-diagonal values in column j.
		for i := 0; i < nrows; i++ {
			if i == j { continue; }
			cij := c.Data[i][j]
			for col := 0; col < ncols; col++ {
				c.Data[i][col] -= cij * c.Data[j][col]
			}
		}
	}
	return c, nil
}

