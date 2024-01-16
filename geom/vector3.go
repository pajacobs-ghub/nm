// vector3.go
// A package for working with geometric vectors in three dimensions.
// PJ 2022-04-26

package geom

import (
	"fmt"
	"math"
)

type Vector3 struct {
	X, Y, Z float64
}

func (v Vector3) String() string {
	return fmt.Sprintf("(%0.6f, %0.6f, %0.6f)", v.X, v.Y, v.Z)
}

func (v Vector3) ApproxEquals(other Vector3, tol float64) bool {
	return (math.Abs(v.X-other.X) <= tol) &&
		(math.Abs(v.Y-other.Y) <= tol) &&
		(math.Abs(v.Z-other.Z) <= tol)
}

func (v Vector3) Add(other Vector3) Vector3 {
	return Vector3{v.X + other.X, v.Y + other.Y, v.Z + other.Z}
}

func (v Vector3) Sub(other Vector3) Vector3 {
	return Vector3{v.X - other.X, v.Y - other.Y, v.Z - other.Z}
}

func (v Vector3) Mul(m float64) Vector3 {
	return Vector3{m * v.X, m * v.Y, m * v.Z}
}

func (v Vector3) Dot(other Vector3) float64 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

// Alternative is to pass pointer values.
func (v *Vector3) Dotp(other *Vector3) float64 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}
