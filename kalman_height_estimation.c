#include "stdlib.h"
#include "kalman_height_estimation.h"
#include "CLinearAlgebra.h"

#define STD_DEV_V_H      0.001f // process noise
#define STD_DEV_W_H      0.01f // sensor noise

typedef struct kalman_height_estimation {
	float dt;
	CLinearAlgebra_Matrix_t A;
	CLinearAlgebra_Matrix_t B;
	CLinearAlgebra_Matrix_t C;
	CLinearAlgebra_Matrix_t V;
	CLinearAlgebra_Matrix_t W;
	CLinearAlgebra_Matrix_t X;
	CLinearAlgebra_Matrix_t P;
	CLinearAlgebra_Matrix_t Xprev;
	CLinearAlgebra_Matrix_t Pprev;
	CLinearAlgebra_Matrix_t U;
	CLinearAlgebra_Matrix_t Y;
	CLinearAlgebra_Matrix_t E;
	CLinearAlgebra_Matrix_t S;
	CLinearAlgebra_Matrix_t K;
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
	handle->A.NumRows = 2;
	handle->A.NumCols = 2;
	CLinearAlgebra_Matrix_SetValue(&handle->A, 0, 0, 1.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->A, 0, 1, handle->dt);
	CLinearAlgebra_Matrix_SetValue(&handle->A, 1, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->A, 1, 1, 1.0f);

	handle->B.NumRows = 2;
	handle->B.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->B, 0, 0, 0.5f * handle->dt * handle->dt);
	CLinearAlgebra_Matrix_SetValue(&handle->B, 1, 0, handle->dt);

	handle->C.NumRows = 1;
	handle->C.NumCols = 2;
	CLinearAlgebra_Matrix_SetValue(&handle->C, 0, 0, 1.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->C, 0, 1, 0.0f);

	handle->V.NumRows = 2;
	handle->V.NumCols = 2;
	CLinearAlgebra_Matrix_SetValue(&handle->V, 0, 0, STD_DEV_V_H * STD_DEV_V_H);
	CLinearAlgebra_Matrix_SetValue(&handle->V, 0, 1, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->V, 1, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->V, 1, 1, STD_DEV_V_H * STD_DEV_V_H);

	handle->W.NumRows = 1;
	handle->W.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->W, 0, 0, STD_DEV_W_H * STD_DEV_W_H);

	handle->X.NumRows = 2;
	handle->X.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->X, 0, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->X, 1, 0, 0.0f);

	handle->P.NumRows = 2;
	handle->P.NumCols = 2;
	CLinearAlgebra_Matrix_SetValue(&handle->P, 0, 0, 1.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->P, 0, 1, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->P, 1, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->P, 1, 1, 1.0f);

	handle->Xprev.NumRows = 2;
	handle->Xprev.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->Xprev, 0, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->Xprev, 1, 0, 0.0f);

	handle->Pprev.NumRows = 2;
	handle->Pprev.NumCols = 2;
	CLinearAlgebra_Matrix_SetValue(&handle->Pprev, 0, 0, 1.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->Pprev, 0, 1, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->Pprev, 1, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->Pprev, 1, 1, 1.0f);

	handle->U.NumRows = 1;
	handle->U.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->U, 0, 0, 0.0f);

	handle->Y.NumRows = 1;
	handle->Y.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->Y, 0, 0, 0.0f);

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

	handle->E.NumRows = 1;
	handle->E.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->E, 0, 0, 0.0f);

	handle->S.NumRows = 1;
	handle->S.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->S, 0, 0, 0.0f);

	handle->K.NumRows = 2;
	handle->K.NumCols = 1;
	CLinearAlgebra_Matrix_SetValue(&handle->K, 0, 0, 0.0f);
	CLinearAlgebra_Matrix_SetValue(&handle->K, 1, 0, 0.0f);

	CLinearAlgebra_Matrix_SetValue(&handle->U, 0, 0, accel_z);
	CLinearAlgebra_Matrix_SetValue(&handle->Y, 0, 0, height);


	// X = A*Xprev + B*U
	handle->X = CLinearAlgebra_Matrix_Add(CLinearAlgebra_Matrix_Multiply(handle->A, handle->Xprev),
	                                      CLinearAlgebra_Matrix_Multiply(handle->B, handle->U));

	// Ppri_h = (A_h * Ppost_h * A_h.t()) + V_h;
	handle->P = CLinearAlgebra_Matrix_Add(CLinearAlgebra_Matrix_Multiply(CLinearAlgebra_Matrix_Multiply(handle->A, handle->Pprev),
	                                      CLinearAlgebra_Matrix_Transpose(handle->A)),
	                                      handle->V);

	// E_h = Y_h - (C_h * Xpri_h);
	handle->E = CLinearAlgebra_Matrix_Subtract(handle->Y,
	            CLinearAlgebra_Matrix_Multiply(handle->C, handle->X));

	// S_h = (C_h * Ppri_h * C_h.t()) + W_h;
	handle->S = CLinearAlgebra_Matrix_Add(CLinearAlgebra_Matrix_Multiply(CLinearAlgebra_Matrix_Multiply(handle->C, handle->P),
	                                      CLinearAlgebra_Matrix_Transpose(handle->C)),
	                                      handle->W);

	// K_h = (Ppri_h * C_h.t()) * S_h.inverse();
	handle->K = CLinearAlgebra_Matrix_Multiply(
	                CLinearAlgebra_Matrix_Multiply(handle->P,
	                        CLinearAlgebra_Matrix_Transpose(handle->C)),
	                CLinearAlgebra_Matrix_Inverse(handle->S));

	// Xpost_h = Xpri_h + (K_h * E_h);
	handle->Xprev = CLinearAlgebra_Matrix_Add(handle->X,
	                CLinearAlgebra_Matrix_Multiply(handle->K, handle->E));

	// Ppost_h = Ppri_h - (K_h * S_h * K_h.t());
	handle->Pprev = CLinearAlgebra_Matrix_Subtract(handle->P, CLinearAlgebra_Matrix_Multiply(
	                    CLinearAlgebra_Matrix_Multiply(handle->K, handle->S),
	                    CLinearAlgebra_Matrix_Transpose(handle->K)
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

	*height = CLinearAlgebra_Matrix_GetValue(handle->X, 0, 0);

	return ERR_CODE_SUCCESS;
}