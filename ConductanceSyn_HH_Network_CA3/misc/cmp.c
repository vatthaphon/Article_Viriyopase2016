/*
 * cmp.c
 *
 *  Created on: Oct 10, 2011
 *      Author: viriyopa
 *
 * CMP Two-value comparison
 *	val = cmp(x, y, tol_eq)
 * Input
 *  x		the first number.
 *  y       the second number.
 *  tol_eq  if the first and second numbers are different less than tol_eq, we say that the two numbers are equal.
 * Output
 *  val     0   : two numbers are the same.
 *          -1  : the first number is less than the second number.
 *          1   : the first number is greater than the second number.
 */

#include <math.h>

int cmp(double x, double y, double tol_eq)
{
	if (fabs(x-y) < tol_eq)
	{
		return 0;
	}
	else if (x < y)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}
