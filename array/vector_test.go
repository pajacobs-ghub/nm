/** vector_test.go
 * Try out some of the Vector functions.
 * PJ 2025-03-06
 */

package array

import (
	"testing"
	_ "fmt"
	"math"
)

func TestVector(t *testing.T) {
	v0 := Vector{}                            // Empty Vector.
	// fmt.Println("v0=", v0.String())
	if !v0.IsEmpty() {
		t.Errorf("Vector v0 appears to be not empty.")
	}
	v1 := Vector{[]float64{1.0, 2.0, 3.0}}    // One way to initialize.
	// fmt.Println("v1=", v1.String())
	if len(v1.Data) != 3 || v1.Data[0] != 1.0 {
		t.Errorf("Vector initialize error v1= %v want= (1.0, 2.0, 3.0)", v1.String())
	}
	v2 := VectorZeros(4)                      // Another way.
	// fmt.Println("v2=", v2.String())
	if len(v2.Data) != 4 || v2.Data[0] != 0.0 {
		t.Errorf("Vector initialize error v2= %v want= (0.0, 0.0, 0.0, 0.0)", v2.String())
	}
	v3 := VectorOnes(4)                       // And another way.
	// fmt.Println("v3=", v3.String())
	if len(v3.Data) != 4 || v3.Sum() != 4.0 { // Exactly-representable numbers used.
		t.Errorf("Vector initialize error v3= %v want= (1.0, 1.0, 1.0, 1.0)", v3.String())
	}
	m := v3.Mean()
	if m != 1.0 { // Exactly-representable numbers used.
		t.Errorf("Vector Mean error m= %v want= 1.0", m)
	}
	m2 := v3.Mag()
	if math.Abs(m2 - 2.0) > 1.0e-9 {
		t.Errorf("Vector Mag error m2= %v want= 2.0", m2)
	}
	v3.Normalize()
	if len(v3.Data) != 4 || math.Abs(v3.Sum() - 2.0) > 1.0e-9 {
		t.Errorf("Vector initialize error v3= %v want= (0.5, 0.5, 0.5, 0.5)", v3.String())
	}
	v4 := VectorCopyArray([]float64{1.1, 2.2, 3.3, 4.4})
	v5 := v4.Clone()
	// fmt.Println("v4=", v4.String())
	if len(v4.Data) != 4 || !v4.ApproxEquals(v5, 1.0e-9) {
		t.Errorf("Vector clone error v5= %v want= %v", v5.String(), v4.String())
	}
	v6 := VectorZeros(4)
	v6.Add(&v4, &v5)
	v5.Scale(2.0)
	if len(v6.Data) != 4 || !v6.ApproxEquals(v5, 1.0e-9) {
		t.Errorf("Vector add error v6= %v want= %v", v6.String(), v5.String())
	}
	v6.SetScalar(0.0)
	if len(v6.Data) != 4 || v6.Data[0] != 0.0 {
		t.Errorf("Vector set-zeros error v6= %v want= (0.0, 0.0, 0.0, 0.0)", v6.String())
	}
	v6.Blend(&v4, &v4, 0.5, 1.5)
	if len(v6.Data) != 4 || !v6.ApproxEquals(v5, 1.0e-9) {
		t.Errorf("Vector add-with-scale error v6= %v want= %v", v6.String(), v5.String())
	}
	v8 := VectorOnes(4)
	s := VectorDot(&v8, &v8)
	if math.Abs(s - 4.0) > 1.0e-9 {
		t.Errorf("Vector dot product error s= %v want= %v", s, 4.0)
	}
}
