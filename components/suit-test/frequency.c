#include <stdio.h>
#include <math.h>
#include <string.h>
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "frequency.h"
static const char* TAG = __FILENAME__;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                          F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int Frequency(const unsigned char *bitArray, int bitArrayLen)
{
	int i;
	double sum, s_obs, f, p_value, sqrt2 = 1.41421356237309504880;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t\t      FREQUENCY TEST\n");
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
#endif

	if ( bitArrayLen < 100 )
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\t   n=%d is too short\n", bitArrayLen);
#endif
		return INPUT_LENGTH_TOO_SHORT;
	}
	
	sum = 0.0;
	for ( i = 0; i < bitArrayLen; i++ )
	{
		sum += 2 * (int)bitArray[i] - 1;
	}
	s_obs = fabs(sum) / sqrt((double)(bitArrayLen));
	f = s_obs / sqrt2;
	p_value = cephes_erfc(f);

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\tCOMPUTATIONAL INFORMATION:\n");
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t(a) The nth partial sum = %d\n", (int)sum);
	ESP_LOGI(TAG, "\t\t(b) S_n/n               = %f\n", sum / bitArrayLen);
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
	ESP_LOGI(TAG, "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
#endif

	if ( p_value < ALPHA )
	{
		return FREQUENCY_TEST_FAIL;
	}
	else
	{
		return 0;
	}
}
