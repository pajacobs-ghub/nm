/** nelmin_test.go
 * Try out the Nelder-Mead simplex optimizer.
 *
 * Peter J, 2024-Jan-15
 */

package nelmin

import (
	"fmt"
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

func obj(x []float64) float64 {
	// Assume 3 parameters.
	z0 := x[0] - 1.0
	z1 := x[1] - 1.0
	z2 := x[2] - 1.0
	return z0*z0 + z1*z1 + z2*z2
}

func TestSimplex(t *testing.T) {
	fmt.Println("Test Simplex functions")
	x0 := []float64{1.0, 2.0, 3.0}
	dx := []float64{0.1, 0.2, 0.3}
	smplx, err := MakeSimplexAboutPoint(obj, x0, dx)
	if err != nil {
		t.Errorf("Failed to make simplex, err: %s", err)
	}
	fmt.Println("smplx=", SimplexToJSON(smplx))
	var vmid *Vertex
	vmid, err = Centroid(smplx, 1)
	fmt.Println("vmid=", vmid)
}

func TestMinimizer(t *testing.T) {
	fmt.Println("Test Minimizer functions")
	m := NewMinimizer(obj)
	x := []float64{1.0, 2.0, 3.0}
	dx := []float64{0.1, 0.2, 0.3}
	err := m.MinimizeFromPoint(x, dx, 100)
	if err != nil {
		t.Errorf("Failed to mimimize from point, err: %s", err)
	}
	fmt.Printf("m=%s\n", m.String())
}
