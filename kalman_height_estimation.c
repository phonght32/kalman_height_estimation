#include "stdlib.h"
#include "kalman_height_estimation.h"
#include "clinalg.h"

#define STD_DEV_V_H      0.05f // process noise
#define STD_DEV_W_H      0.3f // sensor noise

typedef struct kalman_height_estimation {
	float dt;
	clinalg_matrix_t A;
	clinalg_matrix_t B;
	clinalg_matrix_t C;
	clinalg_matrix_t V;
	clinalg_matrix_t W;
	clinalg_matrix_t X;
	clinalg_matrix_t P;
	clinalg_matrix_t Xprev;
	clinalg_matrix_t Pprev;
	clinalg_matrix_t U;
	clinalg_matrix_t Y;
	clinalg_matrix_t E;
	clinalg_matrix_t S;
	clinalg_matrix_t K;
} kalman_height_estimation_t;

kalman_height_estimation_handle_t kalman_height_estimation_init(void)
{
	kalman_height_estimation_handle_t handle = calloc(1, sizeof(kalman_height_estimation_t));

	if (handle == NULL)
	{
		return NULL;
	}

	return handle;
}

err_code_t kalman_height_estimation_set_config(kalman_height_estimation_handle_t handle, kalman_height_estimation_cfg_t config)
{
	handle->dt = config.dt;

	return ERR_CODE_SUCCESS;
}

err_code_t kalman_height_estimation_config(kalman_height_estimation_handle_t handle)
{
	handle->A = clinalg_matrix_create(2, 2, (float[]){ 
		1.0f, handle->dt, 
		0.0f, 1.0f 
	});

	handle->B = clinalg_matrix_create(2, 1, (float[]){
		0.5f * handle->dt * handle->dt,
		handle->dt
	});

	handle->C = clinalg_matrix_create(1, 2, (float[]){
		1.0f, 0.0f
	});

	handle->V = clinalg_matrix_create(2, 2, (float[]){
		STD_DEV_V_H * STD_DEV_V_H,  0.0f,
		0.0f,                       STD_DEV_V_H * STD_DEV_V_H
	});

	handle->W = clinalg_matrix_create(1, 1, (float[]){
		STD_DEV_W_H * STD_DEV_W_H
	});

	handle->X = clinalg_matrix_create(2, 1, (float[]){
		0.0f,
		0.0f
	});

	handle->P = clinalg_matrix_create(2, 2, (float[]){
		1.0f, 0.0f,
		0.0f, 1.0f
	});

	handle->Xprev = clinalg_matrix_create(2, 1, (float[]){
		0.0f,
		0.0f
	});

	handle->Pprev = clinalg_matrix_create(2, 2, (float[]){
		1.0f, 0.0f,
		0.0f, 1.0f
	});

	handle->U = clinalg_matrix_create(1, 1, (float[]){
		0.0f
	});

	handle->Y = clinalg_matrix_create(1, 1, (float[]){
		0.0f
	});

	return ERR_CODE_SUCCESS;
}

err_code_t kalman_height_estimation_update(kalman_height_estimation_handle_t handle, float accel_z, float height)
{
	/* Check if handle structure is NULL */
	if (handle == NULL)
	{
		return ERR_CODE_NULL_PTR;
	}

	handle->E = clinalg_matrix_create(1, 1, (float[]){
		0.0f
	});

	handle->S = clinalg_matrix_create(1, 1, (float[]){
		0.0f
	});

	handle->K = clinalg_matrix_create(2, 1, (float[]){
		0.0f,
		0.0f
	});

	clinalg_matrix_set_member(&handle->U, 0, 0, accel_z);
	clinalg_matrix_set_member(&handle->Y, 0, 0, height);

	// X = A*Xprev + B*U
	handle->X = clinalg_matrix_add(clinalg_matrix_multiply(handle->A, handle->Xprev),
	                                      clinalg_matrix_multiply(handle->B, handle->U));

	// Ppri_h = (A_h * Ppost_h * A_h.t()) + V_h;
	handle->P = clinalg_matrix_add(clinalg_matrix_multiply(clinalg_matrix_multiply(handle->A, handle->Pprev),
	                                      clinalg_matrix_transpose(handle->A)),
	                                      handle->V);

	// E_h = Y_h - (C_h * Xpri_h);
	handle->E = clinalg_matrix_subtract(handle->Y,
	            clinalg_matrix_multiply(handle->C, handle->X));

	// S_h = (C_h * Ppri_h * C_h.t()) + W_h;
	handle->S = clinalg_matrix_add(clinalg_matrix_multiply(clinalg_matrix_multiply(handle->C, handle->P),
	                                      clinalg_matrix_transpose(handle->C)),
	                                      handle->W);

	// K_h = (Ppri_h * C_h.t()) * S_h.inverse();
	handle->K = clinalg_matrix_multiply(
	                clinalg_matrix_multiply(handle->P,
	                        clinalg_matrix_transpose(handle->C)),
	                clinalg_matrix_inverse(handle->S));

	// Xpost_h = Xpri_h + (K_h * E_h);
	handle->Xprev = clinalg_matrix_add(handle->X,
	                clinalg_matrix_multiply(handle->K, handle->E));

	// Ppost_h = Ppri_h - (K_h * S_h * K_h.t());
	handle->Pprev = clinalg_matrix_subtract(handle->P, clinalg_matrix_multiply(
	                    clinalg_matrix_multiply(handle->K, handle->S),
	                    clinalg_matrix_transpose(handle->K)
	                ));

	return ERR_CODE_SUCCESS;
}

err_code_t kalman_height_estimation_get_height(kalman_height_estimation_handle_t handle, float *height)
{
	/* Check if handle structure is NULL */
	if (handle == NULL)
	{
		return ERR_CODE_NULL_PTR;
	}

	*height = clinalg_matrix_get_member(handle->X, 0, 0);

	return ERR_CODE_SUCCESS;
}
