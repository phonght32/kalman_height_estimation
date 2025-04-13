// MIT License

// Copyright (c) 2025 phonght32

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef __KALMAN_HEIGHT_ESTIMATION_H__
#define __KALMAN_HEIGHT_ESTIMATION_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "err_code.h"

typedef struct kalman_height_estimation *kalman_height_estimation_handle_t;

typedef struct {
	float dt;
} kalman_height_estimation_cfg_t;

kalman_height_estimation_handle_t kalman_height_estimation_init(void);
err_code_t kalman_height_estimation_set_config(kalman_height_estimation_handle_t handle, kalman_height_estimation_cfg_t config);
err_code_t kalman_height_estimation_config(kalman_height_estimation_handle_t handle);
err_code_t kalman_height_estimation_update(kalman_height_estimation_handle_t handle, float accel_z, float height);
err_code_t kalman_height_estimation_get_height(kalman_height_estimation_handle_t handle, float *height);

#ifdef __cplusplus
}
#endif

#endif /* __KALMAN_HEIGHT_ESTIMATION_H__ */
