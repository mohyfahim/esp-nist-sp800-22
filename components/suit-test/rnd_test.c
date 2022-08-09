/**************************************************
* File name: rnd_test.c
* Author: HAN Wei
* Author's blog: http://blog.csdn.net/henter/
* License: MIT License (https://mit-license.org/)
           Copyright (c) 2019 HAN Wei
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "errorCode.h"
#include "c_file_operation.h"
#include "basicDefinition.h"
#include "rnd_test.h"
#include "frequency.h"
#include "blockFrequency.h"
#include "runs.h"
#include "longestRunOfOnes.h"
#include "rank.h"
#include "discreteFourierTransform.h"
#include "nonOverlappingTemplateMatchings.h"
#include "overlappingTemplateMatchings.h"
#include "universal.h"
#include "linearComplexity.h"
#include "serial.h"
#include "approximateEntropy.h"
#include "cumulativeSums.h"
#include "randomExcursions.h"
#include "randomExcursionsVariant.h"

static const char *TAG = "rnd_test";
void convertToBitArray(const unsigned char *input,
                       int inputLen,
                       unsigned char *bitArray)
{
    int i;
    const unsigned char *p;
    unsigned char *q;

    p = input;
    q = bitArray;
    for (i = 0; i < inputLen; i++)
    {
        *q = (*p & 0x80) >> 7;
        q++;
        *q = (*p & 0x40) >> 6;
        q++;
        *q = (*p & 0x20) >> 5;
        q++;
        *q = (*p & 0x10) >> 4;
        q++;
        *q = (*p & 0x8) >> 3;
        q++;
        *q = (*p & 0x4) >> 2;
        q++;
        *q = (*p & 0x2) >> 1;
        q++;
        *q = *p & 0x1;
        q++;

        p++;
    }
}

int randomnessTest(unsigned char *input, int inputLen, int errorCode[], FILE* logfd)
{

    int testResult[15], bitStreamLen, passItemCount;
    unsigned char *bitStream;
    int blockFrequencyBlockLen, nonOverlappingTemplateLen, OverlappingTemplateLen,
        linearComplexityBlockLen, serialBlockLen, approximateEntropyBlockLen;
    double passRate;
    enum testItem
    {
        frequencyFlag,
        blockFrequencyFlag,
        runsFlag,
        longestRunOfOnesFlag,
        matrixRankFlag,
        discreteFourierTransformFlag,
        nonOverlappingTemplateMatchingsFlag,
        overlappingTemplateMatchingsFlag,
        universalFlag,
        linearComplexityFlag,
        serialFlag,
        approximateEntropyFlag,
        cumulativeSumsFlag,
        randomExcursionsFlag,
        randomExcursionsVariantFlag
    } rndTestItem;

    blockFrequencyBlockLen = 128;
    nonOverlappingTemplateLen = 10;
    OverlappingTemplateLen = 10;
    linearComplexityBlockLen = 500;
    serialBlockLen = 4;
    approximateEntropyBlockLen = 10;

#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n\t\t\tParameters config\n");
    fprintf(logfd, "\t\t--------------------------------------------\n");
    fprintf(logfd, "\t\t  blockFrequencyBlockLen = %d\n", blockFrequencyBlockLen);
    fprintf(logfd, "\t\t  nonOverlappingTemplateLen = %d\n", nonOverlappingTemplateLen);
    fprintf(logfd, "\t\t  OverlappingTemplateLen = %d\n", OverlappingTemplateLen);
    fprintf(logfd, "\t\t  linearComplexityBlockLen = %d\n", linearComplexityBlockLen);
    fprintf(logfd, "\t\t  serialBlockLen = %d\n", serialBlockLen);
    fprintf(logfd, "\t\t  approximateEntropyBlockLen = %d\n", approximateEntropyBlockLen);
    fprintf(logfd, "\n");
#endif

    bitStreamLen = inputLen * 8;
    if (!(bitStream = (unsigned char *)INTERNAL_MALLOC(bitStreamLen)))
    {
#ifdef INTERNAL_DEBUG
        fprintf(logfd, "Memory allocation failed!\n");
#endif
        return MEMORY_ALLOCATION_FAIL;
    }
    convertToBitArray(input, inputLen, bitStream);

    rndTestItem = frequencyFlag;
    if ((testResult[rndTestItem] = Frequency(bitStream, bitStreamLen)))
    {
        errorCode[rndTestItem] = testResult[rndTestItem];
#ifdef INTERNAL_DEBUG
        ESP_LOGE(TAG, "Frequency test failed!\n");
        ESP_LOGE(TAG, "Test program exits because all subsequent tests depend on the passing of this test!\n");
#endif
        free(bitStream);
        return errorCode[rndTestItem];
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = BlockFrequency(bitStream,
                                                  bitStreamLen,
                                                  blockFrequencyBlockLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = Runs(bitStream, bitStreamLen)))
    {
        
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = LongestRunOfOnes(bitStream, bitStreamLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = Rank(bitStream, bitStreamLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = DiscreteFourierTransform(bitStream, bitStreamLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = NonOverlappingTemplateMatchings(bitStream,
                                                                   bitStreamLen,
                                                                   nonOverlappingTemplateLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = OverlappingTemplateMatchings(bitStream,
                                                                bitStreamLen,
                                                                OverlappingTemplateLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = Universal(bitStream, bitStreamLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = LinearComplexity(bitStream,
                                                    bitStreamLen,
                                                    linearComplexityBlockLen)))
    {
      
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = Serial(bitStream,
                                          bitStreamLen,
                                          serialBlockLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
        
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = ApproximateEntropy(bitStream,
                                                      bitStreamLen,
                                                      approximateEntropyBlockLen)))
    {
         errorCode[rndTestItem] = testResult[rndTestItem];
       
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = CumulativeSums(bitStream, bitStreamLen)))
    {
      
            errorCode[rndTestItem] = testResult[rndTestItem];
      
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = RandomExcursions(bitStream, bitStreamLen)))
    {
       
            errorCode[rndTestItem] = testResult[rndTestItem];
       
    }
    rndTestItem++;
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\n");
#endif

    if ((testResult[rndTestItem] = RandomExcursionsVariant(bitStream, bitStreamLen)))
    {
        
            errorCode[rndTestItem] = testResult[rndTestItem];
       
    }

#ifdef INTERNAL_DEBUG
    rndTestItem = frequencyFlag;
    fprintf(logfd, "\n\n\t\t\tRandomness Test Result:\n");
    fprintf(logfd, "\t-------------------------------------------\n");
    fprintf(logfd, "\t 1. Frequency test:                          %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 2. Block Frequency test:                    %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 3. Runs test:                               %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 4. Longest Run Of Ones test:                %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 5. Matrix Rank test:                        %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 6. Discrete Fourier Transform test:         %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 7. Non-Overlapping Template Matchings test: %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 8. Overlapping Template Matchings test:     %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t 9. Universal test:                          %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t10. Linear Complexity test:                  %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t11. Serial test:                             %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t12. Approximate Entropy test:                %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t13. Cumulative Sums test:                    %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t14. Random Excursions test:                  %s\n", testResult[rndTestItem] ? "Failure" : "Success");
    rndTestItem++;
    fprintf(logfd, "\t15. randomExcursionsVariant test:            %s\n", testResult[rndTestItem] ? "Failure" : "Success");
#endif

    passItemCount = 0;
    for (rndTestItem = frequencyFlag; rndTestItem <= randomExcursionsVariantFlag; rndTestItem++)
    {
        if (!(testResult[rndTestItem]))
        {
            passItemCount++;
        }
    }
    passRate = (double)(passItemCount) / (sizeof(testResult) / sizeof(int));
#ifdef INTERNAL_DEBUG
    fprintf(logfd, "\t---------------------------------------------\n");
    fprintf(logfd, "\t\t\tSummarization:\n");
    fprintf(logfd, "\t\tPass rate = %.4f%%\n", (passRate * 100));
#endif

    free(bitStream);
    return ESP_OK;
}
