#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "runs.h"
static const char* TAG = __FILENAME__;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              R U N S  T E S T 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int Runs(const unsigned char *bitArray, int bitArrayLen)
{
	int S, k;
	double pi, V, erfc_arg, p_value;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t\t\tRUNS TEST\n");
	ESP_LOGI(TAG, "\t\t------------------------------------------\n");
#endif

	S = 0;
	for ( k = 0; k < bitArrayLen; k++ )
	{
		if ( bitArray[k] )
		{
			S++;
		}
	}
	pi = (double)S / (double)bitArrayLen;

	if ( fabs(pi - (double)(0.5)) >= ((double)(2.0) / sqrt((double)bitArrayLen)) )
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
#endif
		return RUNS_TEST_FAIL;
	}
	else
	{
		V = 1;
		for ( k = 1; k < bitArrayLen; k++ )
		{
			if ( bitArray[k] != bitArray[k-1] )
			{
				V++;
			}
		}
	
		erfc_arg = fabs(V - 2.0 * bitArrayLen * pi * (1.0 - pi)) / (2.0 * pi * (1.0 - pi) * sqrt(2.0 * (double)(bitArrayLen)));
		p_value = cephes_erfc(erfc_arg);
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tCOMPUTATIONAL INFORMATION:\n");
		ESP_LOGI(TAG, "\t\t------------------------------------------\n");
		ESP_LOGI(TAG, "\t\t(a) Pi                        = %f\n", pi);
		ESP_LOGI(TAG, "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
		ESP_LOGI(TAG, "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
		ESP_LOGI(TAG, "\t\t    -----------------------   = %f\n", erfc_arg);
		ESP_LOGI(TAG, "\t\t      2 sqrt(2n) pi (1-pi)\n");
		ESP_LOGI(TAG, "\t\t------------------------------------------\n");
#endif
		if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		{
#ifdef INTERNAL_DEBUG
			ESP_LOGI(TAG, "\n\t\tp_value = %f\n", p_value);
			ESP_LOGI(TAG, "WARNING:  P_VALUE IS OUT OF RANGE.\n");
#endif
			return P_VALUE_OUT_OF_RANGE;
		}
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
#endif
	}

	if ( p_value < ALPHA )
	{
		return RUNS_TEST_FAIL;
	}
	else
	{
		return 0;
	}
}
