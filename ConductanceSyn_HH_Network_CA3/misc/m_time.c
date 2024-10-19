/*
 * m_time.c
 *
 *  Created on: Oct 15, 2011
 *      Author: attha
 */

#include <time.h>
#include <string.h>

char *getCurrentTime()
{
	time_t now;
	time(&now);
	char *tmp = ctime(&now);

	tmp[strlen(tmp)-1]='\0';

	return tmp;
}
