/*
 * Matrix.c
 *
 *  Created on: Feb 6, 2012
 *      Author: viriyopa
 */

/*
 * 0<=x<=(dim(array,1)=:x_l)-1, 0<=y<=(dim(array,2)=:y_l)-1
 */
double getAEle(double array[], int x, int y, int x_l, int y_l)
{
	return array[y+x*y_l];
}

double getIntPtrEle(int *int_ptr, int x, int y, int x_l, int y_l)
{
	return *(int_ptr + y + x*y_l);
}

double getDoublePtrEle(double *double_ptr, int x, int y, int x_l, int y_l)
{
	return *(double_ptr + y+ x*y_l);
}

/*
 * 0<=x<=(dim(array,1)=:x_l)-1, 0<=y<=(dim(array,2)=:y_l)-1
 */
void setAEle(double array[], double ele, int x, int y, int x_l, int y_l)
{
	array[y+x*y_l]=ele;
}

void setIntPtrEle(int *int_ptr, double ele, int x, int y, int x_l, int y_l)
{
	*(int_ptr + y + x*y_l) = ele;
}

void setDoublePtrEle(double *double_ptr, double ele, int x, int y, int x_l, int y_l)
{
	*(double_ptr + y + x*y_l) = ele;
}

