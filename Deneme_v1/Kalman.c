#include <stdio.h>


// x_k,y_k are the measured position in x and y directions at time k.


main()
{
	// Define constant matrices
	int A[4][4];
	int A_transpose[4][4];
	// Define variable matrices
	int X_pred[4][1]; // Predicted vector
	int X[4][1]; // State vector
	double P_pred[4][4]; // Predicted covariance
	double P[4][4]; // Corrected covariance
	// Define constants for the initial prediction and covariance matrices
	int x_center = 0;
	int y_center = 0;

	// Initialize A matrix
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			A[i][j] = (i == j ? 1 : 0);
			A_transpose[i][j] = (i == j ? 1 : 0);
		}
	}

	A[0][2] = 1;
	A[1][3] = 1;
	A_transpose[2][0] = 1;
	A_transpose[3][1] = 1;

	//Initialize P matrix
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			P[i][j] = (i == j ? 100 : 0);
		}
	}

	//Calculate Error covarience prediction (equation 2 in pdf)
	double A_cross_P[4][4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4;j++)
		{
			A_cross_P[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				A_cross_P[i][j] += A[i][k] * P[k][j];
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4;j++)
		{
			P_pred[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				P_pred[i][j] += A_cross_P[i][k] * A_transpose[k][j];
			}
		}
	}
	// Add Q matrix 
	P_pred[0][0] = P_pred[0][0] + 0.01;
	P_pred[1][1] = P_pred[1][1] + 0.01;
	P_pred[2][2] = P_pred[2][2] + 0.01;
	P_pred[3][3] = P_pred[3][3] + 0.01;

	//Calculate Kalman Gain (equation 3 in pdf)
	//Define constant matrices for Kalman gain
	int H[2][4];
	int H_transpose[4][2];
	double R[2][2];

	H[0][0] = 1; H[0][1] = 0; H[0][2] = 1; H[0][3] = 0;
	H[1][0] = 0; H[1][1] = 1; H[1][2] = 0; H[1][3] = 1;

	H_transpose[0][0] = 1;  H_transpose[0][1] = 0;
	H_transpose[1][0] = 0;  H_transpose[1][1] = 1;
	H_transpose[2][0] = 1;  H_transpose[2][1] = 0;
	H_transpose[3][0] = 0;  H_transpose[3][1] = 1;

	R[0][0] = 0.2845; R[0][1] = 0.0045;
	R[1][0] = 0.0045; R[2][0] = 0.0455;
	// Multiply H with Ppred
	double H_cross_Ppred[2][4];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4;j++)
		{
			H_cross_Ppred[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				H_cross_Ppred[i][j] += H[i][k] * P_pred[k][j];
			}
		}
	}

	// Multiply H_cross_Ppred with H transpose
	double result[2][2];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2;j++)
		{
			result[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				result[i][j] += H_cross_Ppred[i][k] * H_transpose[k][j];
			}
		}
	}

	result[0][0] += R[0][0]; result[0][1] += R[0][1];
	result[1][0] += R[1][0]; result[1][1] += R[1][1];
	//Take the inverse of the result
	double determinant = (result[0][0] * result[1][1]) - (result[0][1] * result[1][0]);
	double coeff = 1 / determinant;
	result[0][0] = coeff * result[1][1]; result[0][1] = -(coeff*result[0][1]);
	result[1][0] = -(coeff*result[1][0]); result[1][1] = coeff * result[0][0];

	// Multiply Ppred with H transpose
	double Ppred_cross_H_tranpose[4][2];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2;j++)
		{
			Ppred_cross_H_tranpose[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				Ppred_cross_H_tranpose[i][j] += P_pred[i][k] * H_transpose[k][j];
			}
		}
	}

	//Final multiplication for the Kalman gain Ppred_cross_H_transpose cross result
	double Kalman_gain[4][2];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2;j++)
		{
			Kalman_gain[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				Kalman_gain[i][j] += Ppred_cross_H_tranpose[i][k] * result[k][j];
			}
		}
	}


	// Update State (equation 4 in pdf)
	int Z[2][1];
	Z[0][0] = 0; // measured x position
	Z[1][0] = 0; // measured y position

	// H_cross_Xpred
	int H_cross_Xpred[2][1];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 1;j++)
		{
			H_cross_Xpred[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				H_cross_Xpred[i][j] += H[i][k] * X_pred[k][j];
			}
		}
	}

	int Z_HXpred[2][1];
	Z_HXpred[0][0] = Z[0][0] - H_cross_Xpred[0][0];
	Z_HXpred[1][0] = Z[1][0] - H_cross_Xpred[1][0];
	// Kalman gain cross Z_HXpred
	double K_cross_Z_HXpred[4][1];


	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 1;j++)
		{
			K_cross_Z_HXpred[i][j] = 0;
			for (int k = 0; k < 2; k++)
			{
				K_cross_Z_HXpred[i][j] += Kalman_gain[i][k] * Z_HXpred[k][j];
			}
		}
	}

	X[0][0] = X_pred[0][0] + K_cross_Z_HXpred[0][0];
	X[1][0] = X_pred[1][0] + K_cross_Z_HXpred[1][0];
	X[2][0] = X_pred[2][0] + K_cross_Z_HXpred[2][0];
	X[3][0] = X_pred[3][0] + K_cross_Z_HXpred[3][0];

	// Calculate the Error covariance (equation 5 in pdf) (I-KH)Ppred
	//K cross H
	double K_cross_H[4][4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4;j++)
		{
			K_cross_H[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				K_cross_H[i][j] += Kalman_gain[i][k] * H[k][j];
			}
		}
	}
	// I - KH
	double result2[4][4];
	result2[0][0] = 1 - K_cross_H[0][0]; result2[0][1] =  -K_cross_H[0][1]; result2[0][2] = - K_cross_H[0][0]; result2[0][3] = - K_cross_H[0][0];
	result2[1][0] = - K_cross_H[1][0]; result2[1][1] = 1 - K_cross_H[1][1]; result2[1][2] = - K_cross_H[0][0]; result2[1][3] = - K_cross_H[0][0];
	result2[2][0] = - K_cross_H[2][0]; result2[0][0] = - K_cross_H[2][1]; result2[2][2] = 1 - K_cross_H[0][0]; result2[2][3] = - K_cross_H[0][0];
	result2[3][0] = - K_cross_H[3][0]; result2[0][0] = - K_cross_H[3][1]; result2[3][2] = - K_cross_H[0][0]; result2[3][3] = 1 - K_cross_H[0][0];
	// (I-kH)* Ppred
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4;j++)
		{
			P[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				P[i][j] += result2[i][k] * P_pred[k][j];
			}
		}
	}

}
