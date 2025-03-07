// vector.go
// A package for working with vectors of numbers of arbitrary length.
// The intention is to make the arithmetic code in nelmin a bit cleaner.
//
// PJ 2025-03-06

package array

import (
	"fmt"
	"math"
	"bytes"
	"errors"
)

type Vector struct {
	Data []float64
}

func (v Vector) IsEmpty() bool {
	return len(v.Data) == 0
}

func NewVector(n int) *Vector {
	return &Vector{Data: make([]float64, n)}
}

func NewVectorFromArray(data []float64) *Vector {
	n := len(data)
	z := Vector{Data: make([]float64, n)}
	for i := 0; i < n; i++ {
		z.Data[i] = data[i]
	}
	return &z
}

func (a *Vector) Clone() *Vector {
	n := len(a.Data)
	z := Vector{Data: make([]float64, n)}
	for i := 0; i < n; i++ {
		z.Data[i] = a.Data[i]
	}
	return &z
}

// Generate a JSON-compatible string representation.
func (a *Vector) String() string {
	var b bytes.Buffer
	n := len(a.Data)
	b.WriteString("[")
	for i, d := range a.Data {
		b.WriteString(fmt.Sprintf("%g", d))
		if i+1 < n {
			b.WriteString(", ")
		}
	}
	b.WriteString("]")
	return b.String()
}

func (z *Vector) SetFromScalar(a float64) *Vector {
	n := len(z.Data)
	if n == 0 {
		return z
	}
	for i := 0; i < n; i++ {
		z.Data[i] = a
	}
	return z
}

func (z *Vector) SetFromVector(a Vector) (*Vector, error) {
	n := len(z.Data)
	if n != len(a.Data) {
		msg := fmt.Sprintf("Inconsistent array lengths z:%v a:%v", len(z.Data), len(a.Data))
		return z, errors.New(msg)
	}
	if n == 0 {
		return z, nil
	}
	for i := 0; i < n; i++ {
		z.Data[i] = a.Data[i]
	}
	return z, nil
}

func (a *Vector) Sum() float64 {
	s := 0.0
	for _, d := range a.Data {
		s += d
	}
	return s
}

func (a *Vector) Mean() float64 {
	n := len(a.Data)
	if n == 0 {
		return 0.0
	}
	s := 0.0
	for _, d := range a.Data {
		s += d
	}
	return s/float64(n)
}

// Euclidian (L2) norm
func (a *Vector) Mag() float64 {
	n := len(a.Data)
	if n == 0 {
		return 0.0
	}
	s := 0.0
	for i := 0; i < n; i++ {
		s += a.Data[i] * a.Data[i]
	}
	return math.Sqrt(s)
}

func (z *Vector) Normalize() *Vector {
	mag := z.Mag()
	if mag == 0.0 {
		return z
	}
	n := len(z.Data)
	for i := 0; i < n; i++ {
		z.Data[i] /= mag
	}
	return z
}

func (a *Vector) ApproxEquals(other *Vector, tol float64) bool {
	n := len(a.Data)
	if n != len(other.Data) {
		return false
	}
	for i := 0; i < n; i++ {
		aa := a.Data[i]
		bb := other.Data[i]
	    // Relative comparison for large numbers, absolute comparison for small numbers
		if math.Abs(aa-bb)/(0.5*(math.Abs(aa)+math.Abs(bb)+1.0)) > tol {
			return false
		}
	}
	return true
}

func (z *Vector) Scale(a float64) *Vector {
	for i := 0; i < len(z.Data); i++ {
		z.Data[i] *= a
	}
	return z
}

// For the arithmetic function signatures, use the math/big package as a model.
// If results are always pre-allocated, we should have better control
// of the memory required for our expressions.
// Also, we allow aliasing of the arguments so that we can achieve certain
// effects, e.g. z = z + a can be obtained as z = z.Add(z,a)

func (z *Vector) Add(a, b *Vector) (*Vector, error) {
	n := len(z.Data)
	if n != len(a.Data) || n != len(b.Data) {
		msg := fmt.Sprintf("Inconsistent array lengths z:%v a:%v b:%v",
			len(z.Data), len(a.Data), len(b.Data))
		return z, errors.New(msg)
	}
	for i := 0; i < n; i++ {
		z.Data[i] = a.Data[i] + b.Data[i]
	}
	return z, nil
}

func (z *Vector) Blend(a *Vector, b *Vector, sa float64, sb float64) (*Vector, error) {
	n := len(z.Data)
	if n != len(a.Data) || n != len(b.Data) {
		msg := fmt.Sprintf("Inconsistent array lengths z:%v a:%v b:%v",
			len(z.Data), len(a.Data), len(b.Data))
		return z, errors.New(msg)
	}
	for i := 0; i < n; i++ {
		z.Data[i] = sa*a.Data[i] + sb*b.Data[i]
	}
	return z, nil
}

func (z *Vector) Sub(a, b *Vector) (*Vector, error) {
	n := len(z.Data)
	if n != len(a.Data) || n != len(b.Data) {
		msg := fmt.Sprintf("Inconsistent array lengths z:%v a:%v b:%v",
			len(z.Data), len(a.Data), len(b.Data))
		return z, errors.New(msg)
	}
	for i := 0; i < n; i++ {
		z.Data[i] = a.Data[i] - b.Data[i]
	}
	return z, nil
}

func VectorDot(a, b *Vector) (float64, error) {
	n := len(a.Data)
	if n != len(b.Data) {
		msg := fmt.Sprintf("Inconsistent array lengths a:%v b:%v",
			len(a.Data), len(b.Data))
		return 0.0, errors.New(msg)
	}
	s := 0.0
	for i := 0; i < n; i++ {
		s += a.Data[i] * b.Data[i]
	}
	return s, nil
}
