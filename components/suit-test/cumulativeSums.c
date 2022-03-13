#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "cumulativeSums.h"
static const char* TAG = __FILENAME__;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		    C U M U L A T I V E  S U M S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int CumulativeSums(const unsigned char *bitArray,
                   int bitArrayLen)
{
	int S, sup, inf, z = 0, zrev = 0, k;
	double sum1, sum2, p_value1, p_value2;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
	ESP_LOGI(TAG, "\t\t-------------------------------------------\n");
#endif

	if ( bitArrayLen < 100)
	{
#ifdef INTERNAL_DEBUG
        ESP_LOGI(TAG, "\t\tError: The minimal input length in bits is 100!\n");
		ESP_LOGI(TAG, "\t\t--------------------------------------------\n");
#endif
        return INPUT_LENGTH_TOO_SHORT;      
	}

	S = 0;
	sup = 0;
	inf = 0;
	for ( k = 0; k < bitArrayLen; k++ )
    {
		bitArray[k] ? S++ : S--;
		if ( S > sup )
        {
            sup++;    /* wave peak */
        }
		if ( S < inf )
        {
            inf--;    /* wave trough */
        }
		z = (sup > (-inf)) ? sup : (-inf);    /* the largest absolute values of set: { sup, inf } */
		zrev = ((sup - S) > (S - inf)) ? (sup - S) : (S - inf);
	}
	
	// forward
	sum1 = 0.0;
	for ( k = ((-bitArrayLen) / z + 1) / 4; k <= (bitArrayLen / z - 1) / 4; k++ )
    {
		sum1 += cephes_normal(((4 * k + 1) * z) / sqrt((double)(bitArrayLen)));
		sum1 -= cephes_normal(((4 * k - 1) * z) / sqrt((double)(bitArrayLen)));
	}

	sum2 = 0.0;
	for ( k = ((-bitArrayLen) / z - 3) / 4; k <= (bitArrayLen / z - 1) / 4; k++ )
    {
		sum2 += cephes_normal(((4 * k + 3) * z) / sqrt((double)(bitArrayLen)));
		sum2 -= cephes_normal(((4 * k + 1) * z) / sqrt((double)(bitArrayLen)));
	}

	p_value1 = 1.0 - sum1 + sum2;
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\tCOMPUTATIONAL INFORMATION:\n");
	ESP_LOGI(TAG, "\t\t-------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t(a) The maximum partial sum = %d\n", z);
	ESP_LOGI(TAG, "\t\t-------------------------------------------\n");
#endif
	if ( isNegative(p_value1) || isGreaterThanOne(p_value1) )
    {
#ifdef INTERNAL_DEBUG
        ESP_LOGI(TAG, "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
#endif
        return P_VALUE_OUT_OF_RANGE;
    }
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "%s\t\tp_value = %f\n\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
#endif

	// backwards
	sum1 = 0.0;
	for ( k = ( (-bitArrayLen) / zrev + 1) / 4; k <= ( bitArrayLen / zrev - 1) / 4; k++ )
    {
		sum1 += cephes_normal(((4 * k + 1) * zrev) / sqrt((double)(bitArrayLen)));
		sum1 -= cephes_normal(((4 * k - 1) * zrev) / sqrt((double)(bitArrayLen)));
	}
	sum2 = 0.0;
	for ( k = ( (-bitArrayLen) / zrev - 3) / 4; k <= (bitArrayLen / zrev - 1) / 4; k++ )
    {
		sum2 += cephes_normal(((4 * k + 3) * zrev) / sqrt((double)(bitArrayLen)));
		sum2 -= cephes_normal(((4 * k + 1) * zrev) / sqrt((double)(bitArrayLen)));
	}
	p_value2 = 1.0 - sum1 + sum2;
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
	ESP_LOGI(TAG, "\t\t-------------------------------------------\n");
	ESP_LOGI(TAG, "\t\tCOMPUTATIONAL INFORMATION:\n");
	ESP_LOGI(TAG, "\t\t-------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t(a) The maximum partial sum = %d\n", zrev);
	ESP_LOGI(TAG, "\t\t-------------------------------------------\n");
#endif
	if ( isNegative(p_value2) || isGreaterThanOne(p_value2) )
    {
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
#endif
        return P_VALUE_OUT_OF_RANGE;
    }
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "%s\t\tp_value = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
#endif

    if ( (p_value1 < ALPHA) || (p_value2 < ALPHA) )
	{
		return CUMULATIVE_SUMS_TEST_FAIL;
	}
	else
	{
		return 0;
	}
}
