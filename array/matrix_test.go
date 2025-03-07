// matrix_test.go
// Try out the simple Matrix functions.
// Peter J. 2025-03-07
//

package array

import (
	"testing"
	"fmt"
	"math"
)

func TestMatrix(t *testing.T) {
	m0 := Matrix{}  // Empty Matrix.
	fmt.Println("m0=", m0.String())
	if !m0.IsEmpty() {
		t.Errorf("Matrix m0 appears to be not empty.")
	}
	m1 := NewMatrix(2,2)
	fmt.Println("m1=", m1.String())
	m2 := NewMatrixFromArray([][]float64{{1.0, 2.0}, {3.0, 4.0}})
	fmt.Println("m2=", m2.String())
	normi := m2.NormInf()
	if math.Abs(normi - 7.0) > 1.0e-9 {
		t.Errorf("Incorrect infinity norm for m2 got=%g want=7.0", normi)
	}
	m3 := NewMatrixFromArray([][]float64{{0.0,2.0,1.0,0.0},{2.0,2.0,0.0,1.0}})
	m4, err := m3.GaussJordanElimination()
	fmt.Println("m4=", m4.String(), " err=", err)
}
