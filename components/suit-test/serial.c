#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "serial.h"
static const char* TAG = __FILENAME__;

/* The function psi2( ) in NIST SP 800-22 test suit may be wrong.
   I wrote a similar function here. */
int psi2(const unsigned char *bitArray,
		 int bitArrayLen,
		 int blockLen,
		 double *sum)
{
	int powLen, i, j, index;
	unsigned int *counterArray;

	powLen = (int)pow((double)(2.0), blockLen);
	if ((counterArray = (unsigned int *)calloc(powLen, sizeof(unsigned int))) == NULL)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "Insufficient memory available.\n");
#endif
		return MEMORY_ALLOCATION_FAIL;
	}

	for (i = 0; i < powLen; i++)
	{
		counterArray[i] = 0;
	}

	for (i = 0; i < bitArrayLen; i++)
	{
		index = bitArray[i];
		for (j = 1; j < blockLen; j++)
		{
			if (bitArray[(i + j) % bitArrayLen] == 0)
			{
				index <<= 1;
			}
			else
			{
				index = (index << 1) | 1;
			}
		}
		counterArray[index]++;
	}

	*sum = 0.0;
	for (i = 0; i < powLen; i++)
	{
		*sum += pow((double)(counterArray[i]), (double)(2));
	}
	*sum = (*sum) * pow((double)(2.0), (double)(blockLen)) / (double)bitArrayLen - (double)bitArrayLen;

	free(counterArray);
	return 0;
}

int Serial(const unsigned char *bitArray,
		   int bitArrayLen,
		   int blockLen)
{
	int errorCode, blockLenUpperBound;
	double psim0, psim1, psim2, del1, del2, p_value1, p_value2;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t       SERIAL TEST\n");
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t COMPUTATIONAL INFORMATION:		  \n");
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t(a) Block length    (m) = %d\n", blockLen);
	ESP_LOGI(TAG, "\t\t(b) Sequence length (n) = %d\n", bitArrayLen);
#endif

	if (blockLen <= 3)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
		ESP_LOGI(TAG, "\t\tWARNING: m <= 3, TEST NOT APPLICABLE.\n");
		ESP_LOGI(TAG, "\t\t\t  m is too small.\n");
#endif
		return INPUT_INVALID_PARAMETER;
	}

	blockLenUpperBound = (int)(floor(log((double)(bitArrayLen)) / log((double)(2.0))) - 2);
	if (blockLen >= blockLenUpperBound)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tNote: The blockSize = %d exceeds recommended value of %d\n",
				blockLen, MAX(1, blockLenUpperBound));
		ESP_LOGI(TAG, "\t\tResults are inaccurate!\n");
		ESP_LOGI(TAG, "\t\t--------------------------------------------\n");
#endif
		return INPUT_INVALID_PARAMETER;
	}

	if ((errorCode = psi2(bitArray, bitArrayLen, blockLen, &psim0)))
	{
		return errorCode;
	}
	if ((errorCode = psi2(bitArray, bitArrayLen, (blockLen - 1), &psim1)))
	{
		return errorCode;
	}
	if ((errorCode = psi2(bitArray, bitArrayLen, (blockLen - 2), &psim2)))
	{
		return errorCode;
	}

	del1 = psim0 - psim1;
	del2 = psim0 - 2.0 * psim1 + psim2;
	p_value1 = cephes_igamc(pow((double)(2.0), (double)(blockLen - 2)), del1 / (double)(2.0));
	p_value2 = cephes_igamc(pow((double)(2.0), (double)(blockLen - 3)), del2 / (double)(2.0));
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t(c) Psi_m               = %f\n", psim0);
	ESP_LOGI(TAG, "\t\t(d) Psi_m-1             = %f\n", psim1);
	ESP_LOGI(TAG, "\t\t(e) Psi_m-2             = %f\n", psim2);
	ESP_LOGI(TAG, "\t\t(f) Del_1               = %f\n", del1);
	ESP_LOGI(TAG, "\t\t(g) Del_2               = %f\n", del2);
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");

	ESP_LOGI(TAG, "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
	ESP_LOGI(TAG, "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
#endif

	if ((p_value1 < ALPHA) || (p_value2 < ALPHA))
	{
		return SERIAL_TEST_FAIL;
	}
	else
	{
		return 0;
	}
}
