/** rkf45.go
 *
 * The Runge-Kutta-Fehlberg ODE stepper.
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
)

type WorkSpace struct {
	arrays [7][]float64
}

func NewWorkSpace(n int) *WorkSpace {
	var ws WorkSpace
	for i := 0; i < 7; i++ {
		ws.arrays[i] = make([]float64, n)
	}
	return &ws
}

/**
 * Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
 *
 * Params:
 *     f is a callable function that returns the derivative of y wrt t
 *        The signature of this function is f(t, y, dydt) where
 *        t is a float value, y is an array of number values
 *        and dydt is the array to hold the computed derivatives.
 *     t0: is the starting value of the independent variable
 *     h: the requested step size
 *     y0: an array of starting values for the dependent variables
 *         It is assumed that the y-elements are indexed 0 .. n-1
 *         where n = y0.length
 *     y1: an array of final values of the dependent variables
 *     err: estimates of the errors in the values of y1
 *
 * Returns:
 *     the final value of the dependent variable
 */
func Step(
	f func(float64, []float64, []float64),
	t0 float64, h float64,
	y0 []float64, y1 []float64, err []float64,
	ws *WorkSpace) float64 {
	n := len(y0)
	// Assuming a system of equations, we need arrays for the intermediate data.
	// We also assume that the workspace arrays are of the correct length.
	k1 := ws.arrays[1]
	k2 := ws.arrays[2]
	k3 := ws.arrays[3]
	k4 := ws.arrays[4]
	k5 := ws.arrays[5]
	k6 := ws.arrays[6]
	ytmp := ws.arrays[0]
	// Build up the sample point information as per the text book descriptions.
	// We assign the result of intermediate array expressions to ytmp
	// because that's needed for D.
	f(t0, y0, k1)
	for j := 0; j < n; j++ {
		ytmp[j] = y0[j] + 0.25*h*k1[j]
	}
	f(t0+h/4.0, ytmp, k2)
	for j := 0; j < n; j++ {
		ytmp[j] = y0[j] + 3.0*h*k1[j]/32.0 + 9.0*h*k2[j]/32.0
	}
	f(t0+3.0*h/8.0, ytmp, k3)
	for j := 0; j < n; j++ {
		ytmp[j] = y0[j] + 1932.0*h*k1[j]/2197.0 - 7200.0*h*k2[j]/2197.0 +
			7296.0*h*k3[j]/2197.0
	}
	f(t0+12.0*h/13.0, ytmp, k4)
	for j := 0; j < n; j++ {
		ytmp[j] = y0[j] + 439.0*h*k1[j]/216.0 - 8.0*h*k2[j] +
			3680.0*h*k3[j]/513.0 - 845.0*h*k4[j]/4104.0
	}
	f(t0+h, ytmp, k5)
	for j := 0; j < n; j++ {
		ytmp[j] = y0[j] - 8.0*h*k1[j]/27.0 + 2.0*h*k2[j] -
			3544.0*h*k3[j]/2565.0 + 1859.0*h*k4[j]/4104.0 - 11.0*h*k5[j]/40.0
	}
	f(t0+h/2.0, ytmp, k6)
	// Now, do the integration as a weighting of the sampled data.
	for j := 0; j < n; j++ {
		y1[j] = y0[j] + 16.0*h*k1[j]/135.0 + 6656.0*h*k3[j]/12825.0 +
			28561.0*h*k4[j]/56430.0 - 9.0*h*k5[j]/50.0 + 2.0*h*k6[j]/55.0
		err[j] = h*k1[j]/360.0 - 128.0*h*k3[j]/4275.0 - 2197.0*h*k4[j]/75240.0 +
			h*k5[j]/50.0 + 2.0*h*k6[j]/55.0
		err[j] = math.Abs(err[j])
	}
	return t0 + h
} // end Step()
