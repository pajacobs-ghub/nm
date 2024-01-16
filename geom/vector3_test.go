/** vector3_test.go
 * Try out some of the Vector3 functions.
 * PJ 2024-01-15
 */

package geom

import (
	"testing"
)

func TestVector3(t *testing.T) {
	v1 := Vector3{1.0, 2.0, 3.0}          // One way to initialize.
	v2 := Vector3{X: 1.1, Y: 2.2, Z: 3.3} // Another way to initialize.
	v1 = v1.Add(v2).Mul(2.0).Sub(Vector3{1.0, 1.0, 1.0})
	v3 := Vector3{X: 3.2, Y: 7.4, Z: 11.6}
	if !v1.ApproxEquals(v3, 1.0e-9) {
		t.Errorf("Vector3 Add error v1= %v want= %v", v1, v3)
	}
}
