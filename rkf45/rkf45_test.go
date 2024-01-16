/** rkf45_test.go
 *
 * Try out the Runge-Kutta ODE stepper.
 *
 * Author: Peter J.
 * Version: 2014-Jun-15, adapted from the mech2700 class example
 *          2014-Jul-09, preallocate work arrays and pass them in.
 *          2018-May-26, work with double or complex numbers
 *          2018-May-30, accept the type of the dependent variables as a parameter
 *          2022-May-20, Build as a single-source-file program.
 *          2022-May-23, Go version
 *          2024-Jan-15, Make part of a Go package.
 */

package rkf45

import (
	"math"
	"testing"
)

/** Test system 1
 * Third-order system with a simple analytic solution.
 * Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
 * Except for the zero-based indexing, the notation is
 * chosen to match that in the text.
 */
func testSystem1(t float64, x []float64, dxdt []float64) {
	dxdt[0] = -8.0/3.0*x[0] - 4.0/3.0*x[1] + x[2] + 12.0
	dxdt[1] = -17.0/3.0*x[0] - 4.0/3.0*x[1] + x[2] + 29.0
	dxdt[2] = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0
}

func analyticSolution1(t float64) []float64 {
	x := math.Exp(-3.0*t) / 6.0 * (6.0 - 50.0*math.Exp(t) +
		10.0*math.Exp(2.0*t) + 34.0*math.Exp(3.0*t))
	y := math.Exp(-3.0*t) / 6.0 * (12.0 - 125.0*math.Exp(t) +
		40.0*math.Exp(2.0*t) + 73.0*math.Exp(3.0*t))
	z := math.Exp(-3.0*t) / 6.0 * (14.0 - 200.0*math.Exp(t) +
		70.0*math.Exp(2.0*t) + 116.0*math.Exp(3.0*t))
	return []float64{x, y, z}
}

func TestODEStepper(t *testing.T) {
	x1 := make([]float64, 3)
	err := make([]float64, 3)
	work := NewWorkSpace(3)
	//
	t0 := 0.0
	t1 := 0.0
	nstep := 10000
	h := 1.0 / float64(nstep)
	x0 := []float64{0.0, 0.0, 0.0}
	for i := 0; i < nstep; i++ {
		t1 = Step(testSystem1, t0, h, x0, x1, err, work)
		for j := 0; j < 3; j++ {
			x0[j] = x1[j]
		}
		t0 = t1
	}
	exact := analyticSolution1(t1)
	failFlag := false
	for j := 0; j < 3; j++ {
		if math.Abs(x1[j]-exact[j]) > 1.0e-9 {
			failFlag = true
		}
	}
	if failFlag {
		t.Errorf("rkf45 step error got= %v want= %v", x1, exact)
	}
} // end TestODEStepper()
