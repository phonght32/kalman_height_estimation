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
	handle->A.rows = 2;
	handle->A.cols = 2;
	clinalg_matrix_set(&handle->A, 0, 0, 1.0f);
	clinalg_matrix_set(&handle->A, 0, 1, handle->dt);
	clinalg_matrix_set(&handle->A, 1, 0, 0.0f);
	clinalg_matrix_set(&handle->A, 1, 1, 1.0f);

	handle->B.rows = 2;
	handle->B.cols = 1;
	clinalg_matrix_set(&handle->B, 0, 0, 0.5f * handle->dt * handle->dt);
	clinalg_matrix_set(&handle->B, 1, 0, handle->dt);

	handle->C.rows = 1;
	handle->C.cols = 2;
	clinalg_matrix_set(&handle->C, 0, 0, 1.0f);
	clinalg_matrix_set(&handle->C, 0, 1, 0.0f);

	handle->V.rows = 2;
	handle->V.cols = 2;
	clinalg_matrix_set(&handle->V, 0, 0, STD_DEV_V_H * STD_DEV_V_H);
	clinalg_matrix_set(&handle->V, 0, 1, 0.0f);
	clinalg_matrix_set(&handle->V, 1, 0, 0.0f);
	clinalg_matrix_set(&handle->V, 1, 1, STD_DEV_V_H * STD_DEV_V_H);

	handle->W.rows = 1;
	handle->W.cols = 1;
	clinalg_matrix_set(&handle->W, 0, 0, STD_DEV_W_H * STD_DEV_W_H);

	handle->X.rows = 2;
	handle->X.cols = 1;
	clinalg_matrix_set(&handle->X, 0, 0, 0.0f);
	clinalg_matrix_set(&handle->X, 1, 0, 0.0f);

	handle->P.rows = 2;
	handle->P.cols = 2;
	clinalg_matrix_set(&handle->P, 0, 0, 1.0f);
	clinalg_matrix_set(&handle->P, 0, 1, 0.0f);
	clinalg_matrix_set(&handle->P, 1, 0, 0.0f);
	clinalg_matrix_set(&handle->P, 1, 1, 1.0f);

	handle->Xprev.rows = 2;
	handle->Xprev.cols = 1;
	clinalg_matrix_set(&handle->Xprev, 0, 0, 0.0f);
	clinalg_matrix_set(&handle->Xprev, 1, 0, 0.0f);

	handle->Pprev.rows = 2;
	handle->Pprev.cols = 2;
	clinalg_matrix_set(&handle->Pprev, 0, 0, 1.0f);
	clinalg_matrix_set(&handle->Pprev, 0, 1, 0.0f);
	clinalg_matrix_set(&handle->Pprev, 1, 0, 0.0f);
	clinalg_matrix_set(&handle->Pprev, 1, 1, 1.0f);

	handle->U.rows = 1;
	handle->U.cols = 1;
	clinalg_matrix_set(&handle->U, 0, 0, 0.0f);

	handle->Y.rows = 1;
	handle->Y.cols = 1;
	clinalg_matrix_set(&handle->Y, 0, 0, 0.0f);

	return ERR_CODE_SUCCESS;
}

err_code_t kalman_height_estimation_update(kalman_height_estimation_handle_t handle,
        float accel_z, float height)
{
	/* Check if handle structure is NULL */
	if (handle == NULL)
	{
		return ERR_CODE_NULL_PTR;
	}

	handle->E.rows = 1;
	handle->E.cols = 1;
	clinalg_matrix_set(&handle->E, 0, 0, 0.0f);

	handle->S.rows = 1;
	handle->S.cols = 1;
	clinalg_matrix_set(&handle->S, 0, 0, 0.0f);

	handle->K.rows = 2;
	handle->K.cols = 1;
	clinalg_matrix_set(&handle->K, 0, 0, 0.0f);
	clinalg_matrix_set(&handle->K, 1, 0, 0.0f);

	clinalg_matrix_set(&handle->U, 0, 0, accel_z);
	clinalg_matrix_set(&handle->Y, 0, 0, height);


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

	*height = clinalg_matrix_get(handle->X, 0, 0);

	return ERR_CODE_SUCCESS;
}
