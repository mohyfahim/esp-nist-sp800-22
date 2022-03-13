#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "errorCode.h"
#include "basicDefinition.h"
#include "cephes.h"
#include "overlappingTemplateMatchings.h" 
static const char* TAG = __FILENAME__;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double Pr(int u, double eta)
{
	int i;
	double sum, p;
	
	if ( u == 0 )
	{
		p = exp(-eta);
	}
	else
	{
		sum = 0.0;
		for ( i = 1; i <= u; i++ )
		{
			sum += exp(-eta - u * log((double)(2.0)) + i * log(eta) - cephes_lgam(i + 1) + cephes_lgam(u) - cephes_lgam(i)-cephes_lgam(u - i + 1));
		}
		p = sum;
	}
	return p;
}

int OverlappingTemplateMatchings(const unsigned char *bitArray,
                                 int bitArrayLen,
								 int templateLen)
{
	int i, k, W_obs, match;
	double lambda, eta, sum, chi2, p_value;
	int blockLen, blockCount, j, K = 5;
	unsigned int nu[6] = { 0, 0, 0, 0, 0, 0 };
	double pi[6];
	unsigned char *sequence;

#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
	ESP_LOGI(TAG, "\t\t---------------------------------------------\n");
#endif

	if ( bitArrayLen < 1000000 )
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tError: n=%d is too short\n", bitArrayLen);
#endif
		return INPUT_LENGTH_TOO_SHORT;
	}

	if ( (templateLen < 2) || ( templateLen > 21) )
    {
#ifdef INTERNAL_DEBUG
        ESP_LOGI(TAG, "\t\tOverlapping template length = %d\n", templateLen);
        ESP_LOGI(TAG, "\t\tError: The parameter value is invalid!\n");
#endif
        return INPUT_INVALID_PARAMETER;
    }

	blockLen = 1032;
	blockCount = bitArrayLen / blockLen;

	if ( (sequence = (unsigned char *) calloc(templateLen, sizeof(unsigned char))) == NULL )
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "\t\tError: Insufficient memory, overlapping Template Matchings test aborted!\n");
#endif
		return MEMORY_ALLOCATION_FAIL;
	}
	
	for ( i = 0; i < templateLen; i++ )
	{
		sequence[i] = 1;
	}
	
	lambda = (double)(blockLen - templateLen + 1) / pow((double)(2.0), templateLen);
	eta = lambda / 2.0;
	sum = 0.0;
	for ( i = 0; i < K; i++ )
	{	/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	for ( i = 0; i < blockCount; i++ )
	{
		W_obs = 0;
		for ( j = 0; j < (blockLen - templateLen + 1); j++ )
		{
			match = 1;
			for ( k = 0; k < templateLen; k++ )
			{
				if ( sequence[k] != bitArray[i * blockLen + j + k] )
				{
					match = 0;
				}
			}
			if ( match == 1 )
			{
				W_obs++;
			}
		}
		if ( W_obs <= 4 )
			nu[W_obs]++;
		else
			nu[K]++;
	}

	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i = 0; i < (K + 1); i++ )
	{
		chi2 += pow((double)nu[i] - (double)blockCount * pi[i], (double)(2.0)) / ((double)blockCount * pi[i]);
	}
	p_value = cephes_igamc((double)(K) / 2.0, chi2 / 2.0);
#ifdef INTERNAL_DEBUG
	ESP_LOGI(TAG, "\t\tCOMPUTATIONAL INFORMATION:\n");
	ESP_LOGI(TAG, "\t\t-----------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t(a) n (sequence_length)      = %d\n", bitArrayLen);
	ESP_LOGI(TAG, "\t\t(b) m (block length of 1s)   = %d\n", templateLen);
	ESP_LOGI(TAG, "\t\t(c) M (length of substring)  = %d\n", blockLen);
	ESP_LOGI(TAG, "\t\t(d) N (number of substrings) = %d\n", blockCount);
	ESP_LOGI(TAG, "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
	ESP_LOGI(TAG, "\t\t(f) eta                      = %f\n", eta);
	ESP_LOGI(TAG, "\t\t-----------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t   F R E Q U E N C Y\n");
	ESP_LOGI(TAG, "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
	ESP_LOGI(TAG, "\t\t-----------------------------------------------\n");
	ESP_LOGI(TAG, "\t\t%3d %3d %3d %3d %3d %3d  %f ",
		nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);
	ESP_LOGI(TAG, "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS");
#endif

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
	{
#ifdef INTERNAL_DEBUG
		ESP_LOGI(TAG, "WARNING:  P_VALUE IS OUT OF RANGE.\n");
#endif
        free(sequence);
		return P_VALUE_OUT_OF_RANGE;
	}

	free(sequence);
	if ( p_value < ALPHA )
	{
		return OVERLAPPING_TEMPLATE_MATCHING_TEST_FAIL;
	}
	else
	{
		return 0;
	}
}

