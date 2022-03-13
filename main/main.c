/**************************************************
* File name: main.c
* Author: Mohy Fahim
* License: MIT License (https://mit-license.org/)
           Copyright (c) 2022 Mohy Fahim
**************************************************/

#include <stdio.h>
#include "errorCode.h"
#include "c_file_operation.h"
#include "rnd_test.h"
#include "esp_err.h"
#include "esp_log.h"
#include "esp_random.h"
#include "bootloader_random.h"

static const char *TAG = "nist-test";

void app_main()
{
    int err_c, fileLen = 32, errorCode[15] = {0};
    uint8_t octetStream[32] = {0};
    ESP_LOGW(TAG, "data debug");
    for (int i = 0; i < 32; i++)
    {
        printf("%u ", octetStream[i]);
    }
    printf("\n");

    bootloader_random_enable();
    esp_fill_random(octetStream, 32);
    bootloader_random_disable();

    for (int i = 0; i < 32; i++)
    {
        printf("%u ", octetStream[i]);
    }
    printf("\n");
    
    if ((err_c = randomnessTest((unsigned char *)(octetStream), fileLen, errorCode)) != ESP_OK)
    {
        ESP_LOGE(TAG, "error on random test : %d ", err_c);
    }
    ESP_LOGW(TAG, "data debug");
    for (int i = 0; i < 32; i++)
    {
        printf("%u ", octetStream[i]);
    }
    printf("\n");
}
