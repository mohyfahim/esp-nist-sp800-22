#include <stdio.h> 
#include <math.h> 
#include <string.h>
#include <stdlib.h>
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "randomExcursionsVariant.h"
static const char* TAG = __FILENAME__;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
            R A N D O M  E X C U R S I O N S  V A R I A N T  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int RandomExcursionsVariant(const unsigned char *bitArray,
                            int bitArrayLen)
{
	int	i, p, J, x, constraint, counter, *S_k;
	int	stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double p_value;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
	ESP_LOGI(TAG, "\t\t--------------------------------------------\n");
#endif

	if ( bitArrayLen < 1000000)
	{
#ifdef INTERNAL_DEBUG
        ESP_LOGI(TAG, "\t\tError: The minimal input length in bits is 1000000!\n");
#endif
        return INPUT_LENGTH_TOO_SHORT;      
	}
	
	if ( (S_k = (int *)calloc(bitArrayLen, sizeof(int))) == NULL )
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tError: Insufficent memory allocated.\n");
#endif
		return MEMORY_ALLOCATION_FAIL;
	}

	J = 0;
	S_k[0] = 2 * (int)bitArray[0] - 1;
	for ( i = 1; i < bitArrayLen; i++ )
	{
		S_k[i] = S_k[i - 1] + 2 * bitArray[i] - 1;
		if ( S_k[i] == 0 )
		{
			J++;
		}
	}

	if ( S_k[bitArrayLen - 1] != 0 )
	{
		J++;
	}
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\tCOMPUTATIONAL INFORMATION:\n");
	ESP_LOGI(TAG, "\t\t--------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t(a) Number Of Cycles (J) = %d\n", J);
	ESP_LOGI(TAG, "\t\t(b) Sequence Length (n)  = %d\n", bitArrayLen);
	ESP_LOGI(TAG, "\t\t--------------------------------------------\n");
#endif

	constraint = MAX( (int)( 0.005 * ( pow( (double)(bitArrayLen), (double)(0.5) ) ) ), (int)(500) );
	if (J < constraint)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\t(c) Rejection Constraint = %d\n", constraint);
		ESP_LOGI(TAG, "\t\tJ < %d\n", constraint);
		ESP_LOGI(TAG, "\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
		ESP_LOGI(TAG, "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
#endif
		free(S_k);
		return RANDOM_EXCURSIONS_VARIANT_TEST_FAIL;
	}

    for ( p = 0; p < 18; p++ )
	{
		x = stateX[p];
		counter = 0;
		for ( i = 0; i < bitArrayLen; i++ )
		{
			if ( S_k[i] == x )
			{
				counter++;
			}
		}
		p_value = cephes_erfc(fabs((double)(counter - J)) / (sqrt((double)(2.0) * J * (4.0 * fabs((double)(x)) - 2))));

		if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		{
#ifdef INTERNAL_DEBUG
			ESP_LOGI(TAG, "\t\t(b) WARNING: P_VALUE IS OUT OF RANGE.\n");
#endif
            free(S_k);
            return P_VALUE_OUT_OF_RANGE;
		}
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
		ESP_LOGI(TAG, "(x = %2d) Total visits = %4d; p-value = %f\n", x, counter, p_value);
#endif
        if ( p_value < ALPHA )
	    {
			free(S_k);
			return RANDOM_EXCURSIONS_VARIANT_TEST_FAIL;
		}
	}

	free(S_k);
    return 0;
}
