#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "c_file_operation.h"
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "nonOverlappingTemplateMatchings.h"
static const char* TAG = __FILENAME__;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		  N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int NonOverlappingTemplateMatchings(const unsigned char *bitArray,
									int bitArrayLen,
									int templateLen)
{
	int numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
							   2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must
	first be constructed, saved into files and then the corresponding
	number of nonperiodic templates for that file be stored in the m-th
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int bit, W_obs;
	FILE *fp = NULL;
	double sum, chi2, lambda, pi[6], varWj;
	double p_value = 0.0;
	int i, j, jj, k, match, SKIP, blockLen, blockCount, K = 5;
	char directory[100];
	unsigned int *Wj = NULL;
	unsigned char *sequence = NULL;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t  NON-OVERLAPPING TEMPLATES MATCHINGS TEST\n");
	ESP_LOGI(TAG, "-------------------------------------------------------------------------------------\n");
#endif

	if ((templateLen < 2) || (templateLen > 21))
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tNon-overlapping template length = %d\n", templateLen);
		ESP_LOGI(TAG, "\t\tError: The parameter value is invalid!\n");
#endif
		return INPUT_INVALID_PARAMETER;
	}

	blockCount = 8;
	blockLen = bitArrayLen / blockCount;

	if ((Wj = (unsigned int *)calloc(blockCount, sizeof(unsigned int))) == NULL)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\tError: Insufficient memory for required work space.\n");
#endif
		return MEMORY_ALLOCATION_FAIL;
	}

	lambda = (blockLen - templateLen + 1) / pow((double)(2.0), templateLen);
	varWj = blockLen * (1.0 / pow((double)(2.0), templateLen) - (2.0 * templateLen - 1.0) / pow((double)(2.0), 2.0 * templateLen));
	sprintf(directory, "templates/template%d", templateLen);

	if ((isNegative(lambda)) || (isZero(lambda)))
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\tError: Lambda (%f) not being positive!\n", lambda);
#endif
		free(Wj);
		return NON_OVERLAPPING_TEMPLATE_MATCHING_TEST_FAIL;
	}

	if ((fp = fopen(directory, "r")) == NULL)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\tError: Template file <%s> not existing\n", directory);
#endif
		free(Wj);
		return FILE_OPERATION_FAIL;
	}

	if ((sequence = (unsigned char *)calloc(templateLen, sizeof(unsigned char))) == NULL)
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\tError: Insufficient memory for required work space.\n");
#endif
		fclose(fp);
		free(Wj);
		return NON_OVERLAPPING_TEMPLATE_MATCHING_TEST_FAIL;
	}

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t  COMPUTATIONAL INFORMATION\n");
	ESP_LOGI(TAG, "-------------------------------------------------------------------------------------\n");
	ESP_LOGI(TAG, "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, blockLen, blockCount,
			templateLen, bitArrayLen);
	ESP_LOGI(TAG, "-------------------------------------------------------------------------------------\n");
	ESP_LOGI(TAG, "\t\tF R E Q U E N C Y\n");
	ESP_LOGI(TAG, "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
	ESP_LOGI(TAG, "-------------------------------------------------------------------------------------\n");
#endif
	if (numOfTemplates[templateLen] < MAX_NUM_OF_TEMPLATES)
	{
		SKIP = 1;
	}
	else
	{
		SKIP = (int)(numOfTemplates[templateLen] / MAX_NUM_OF_TEMPLATES);
	}
	numOfTemplates[templateLen] = (int)numOfTemplates[templateLen] / SKIP;

	sum = 0.0;
	for (i = 0; i < 2; i++)
	{ /* Compute Probabilities */
		pi[i] = exp((-lambda) + i * log(lambda) - cephes_lgam(i + 1));
		sum += pi[i];
	}
	pi[0] = sum;
	for (i = 2; i <= K; i++)
	{ /* Compute Probabilities */
		pi[i - 1] = exp((-lambda) + i * log(lambda) - cephes_lgam(i + 1));
		sum += pi[i - 1];
	}
	pi[K] = 1 - sum;

	for (jj = 0; jj < MIN(MAX_NUM_OF_TEMPLATES, numOfTemplates[templateLen]); jj++)
	{
		sum = 0;

		for (k = 0; k < templateLen; k++)
		{
			fscanf(fp, "%d", &bit);
			sequence[k] = bit;
#ifdef INTERNAL_DEBUG
			ESP_LOGI(TAG, "%d", sequence[k]);
#endif
		}
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, " ");
#endif

		for (i = 0; i < blockCount; i++)
		{
			W_obs = 0;
			for (j = 0; j < (blockLen - templateLen + 1); j++)
			{
				match = 1;
				for (k = 0; k < templateLen; k++)
				{
					if ((int)sequence[k] != (int)bitArray[i * blockLen + j + k])
					{
						match = 0;
						break;
					}
				}
				if (match == 1)
				{
					W_obs++;
					j += templateLen;
				}
			}
			Wj[i] = W_obs;
		}

		sum = 0;
		chi2 = 0.0; /* Compute Chi Square */
		for (i = 0; i < blockCount; i++)
		{
			if (templateLen == 10)
			{
#ifdef INTERNAL_DEBUG
				ESP_LOGI(TAG, "%3d  ", Wj[i]);
#endif
			}
			else
			{
#ifdef INTERNAL_DEBUG
				ESP_LOGI(TAG, "%4d ", Wj[i]);
#endif
			}
			chi2 += (pow(((double)Wj[i] - lambda), (double)(2.0)) / varWj);
		}
		p_value = cephes_igamc(blockCount / 2.0, chi2 / 2.0);

		if (isNegative(p_value) || isGreaterThanOne(p_value))
		{
#ifdef INTERNAL_DEBUG
			ESP_LOGI(TAG, "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");
#endif
			free(sequence);
			fclose(fp);
			free(Wj);
			return P_VALUE_OUT_OF_RANGE;
		}

		if (SKIP > 1)
		{
			fseek(fp, (long)(SKIP - 1) * 2 * templateLen, SEEK_CUR);
		}
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
		ESP_LOGI(TAG, "\n");
#endif
	}

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\n");
#endif
	free(sequence);
	fclose(fp);
	free(Wj);

	if (p_value < ALPHA)
	{
		return NON_OVERLAPPING_TEMPLATE_MATCHING_TEST_FAIL;
	}
	else
	{
		return 0;
	}
}
