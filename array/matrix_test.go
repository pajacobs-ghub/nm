// matrix_test.go
// Try out the simple Matrix functions.
// Peter J. 2025-03-07
//

package array

import (
	"testing"
	_ "fmt"
	"math"
)

func TestMatrix(t *testing.T) {
	m0 := Matrix{}  // Empty Matrix.
	// fmt.Println("m0=", m0.String())
	if !m0.IsEmpty() {
		t.Errorf("Matrix m0 appears to be not empty.")
	}
	m1, err := NewMatrix(2,2)
	// fmt.Println("m1=", m1.String(), " err=", err)
	if err != nil {
		t.Errorf("Did not construct new matrix correctly, m1=%v", m1)
	}
	
	m2, err := NewMatrixFromArray([][]float64{{1.0, 2.0}, {3.0, 4.0}})
	// fmt.Println("m2=", m2.String(), " err=", err)
	normi := m2.NormInf()
	if math.Abs(normi - 7.0) > 1.0e-9 {
		t.Errorf("Incorrect infinity norm for m2 got=%g want=7.0", normi)
	}
	
	m3, err := NewMatrixFromArray([][]float64{{1.0, 2.0}, {3.0}})
	// fmt.Println("m3=", m3.String(), " err=", err)
	if err == nil {
		t.Errorf("Did not detect ragged rows m3=%v", m3)
	}

	m4, err := NewMatrixFromArray([][]float64{{0.0,2.0,1.0,0.0},{2.0,2.0,0.0,1.0}})
	// fmt.Println("before elimination m4=", m4.String(), " err=", err)
	m4, err = m4.GaussJordanElimination()
	// fmt.Println("after elimination m4=", m4.String(), " err=", err)
	m4ref, err := NewMatrixFromArray([][]float64{{1.0, 0.0, -0.5, 0.5},{0.0, 1.0, 0.5, 0.0}})
	// fmt.Println("m4ref=", m4ref.String(), " err=", err)
	if !m4.ApproxEquals(m4ref, 1.0e-9) {
		t.Errorf("Incorrect elimination, m4=%s", m4.String())
	}
	
	m5, _ := NewMatrixFromArray([][]float64{{0.0,2.0,1.0,0.0},{0.0,2.0,0.0,1.0}})
	m6, err := m5.GaussJordanElimination()
	// fmt.Println("m6=", m6.String(), " err=", err)
	if err == nil {
		t.Errorf("Did not detect singular matrix m6=%v", m6)
	}
}
