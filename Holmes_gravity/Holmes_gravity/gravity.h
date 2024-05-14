#pragma once
namespace gravity
{
	void get_gravity_coefficients(std::string filename)
	{
		for (int i = 0; i < M + 1; i++)
		{
			Order2[i] = new int[M + 1];
			for (int j = 0; j < M + 1; j++) Order2[i][j] = 0;
		}

		int i = 0;
		for (int n = 0; n < M + 1; n++)
		{
			for (int m = 0; m <= n; m++)
			{
				Order2[n][m] = i; // заполняем массив индексов
				i++;
			}
		}

		Cnm[0] = 0.0; Cnm[1] = 0.0; Cnm[2] = 0.0; // n = 0; n = 1;
		Snm[0] = 0.0; Snm[1] = 0.0; Snm[2] = 0.0; // n = 0; n = 1;

		i = 3;
		int nn, mm;
		double cc, ss, errc, errs;
		std::ifstream egm(filename);
		if (egm.is_open())
		{
			while (egm >> nn >> mm >> cc >> ss >> errc >> errs)
			{
				Cnm[i] = cc;
				Snm[i] = ss;
				i++;
				if ((nn == M) && (mm == M)) break;
			}
			egm.close();
		}
	}

	// coefficient for recurrent Pnm
	double a(int n, int m)
	{
		if (n < m) return 0.0;
		return sqrt((double)(4 * n * n - 1) / (double)(n * n - m * m));
	}

	// coefficient for recurrent Pnm
	double b(int n, int m)
	{
		if (n < m) return 0.0;
		return sqrt((double)((2 * n + 1)*(n + m - 1)*(n - m - 1)) / (double)((n - m)*(n + m)*(2 * n - 3)));
	}
	
	// coefficient for recurrent P(1)nm
	double f(int n, int m)
	{
		if (n < m) return 0.0;
		return sqrt((double)((n*n - m * m)*(2 * n + 1)) / (double)(2 * n - 1));
	}

	// doesn't work!!! (and currently not in use)
	void P(int n, int m, double cost, double* Ps)
	{
		double p;
		if (n == m)
		{
			Ps[0] = 0.0;
			Ps[1] = 1.0;
		}
		else
		{
			p = a(n, m) * cost * Ps[1] - b(n, m) * Ps[0];
			Ps[0] = Ps[1];
			Ps[1] = p;
		}
	}

	// create all Pnm for given cost = cos(theta)
	void allPs(double cost)
	{
		for (int i = 0; i < M+1; i++)
		{
			for (int j = 0; j < M+1; j++) Ps[i][j] = 0.0;
		}
		Ps[0][0] = 1.0;

		for (int n = 1; n < M+1; n++)
		{
			for (int m = 0; m <= n; m++)
			{
				if (n == m) Ps[n][m] = 1.0;
				else
				{
					if (n > 1) Ps[n][m] = a(n, m) * cost * Ps[n - 1][m] - b(n, m) * Ps[n - 2][m];
					else Ps[n][m] = a(n, m) * cost * Ps[n - 1][m];
				}
			}
		}
	}

	// create all cos(mlambda) and sin(mlambda) for given lambda and M 
	void allml(double sinl, double cosl)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < M+1; j++) ml[i][j] = 0.0;
		}

		ml[0][0] = 0.0; // sin(0)
		ml[1][0] = 1.0; // cos(0)
		ml[0][1] = sinl; // sin(lambda)
		ml[1][1] = cosl; // cos(lambda)

		for (int i = 0; i < 2; i++)
		{
			for (int j = 2; j < M+1; j++) ml[i][j] = 2.0 * cosl * ml[i][j - 1] - ml[i][j - 2];
		}
	}

	/**
		Holmes algorithm of computing gravity acceleration using modified forward column method

		@param r - geocentral distance to the satellite
		@param sint = sin(theta), theta - zenit distance
		@param cost = cos(theta)
		@param sinl = sin(lambda), lambda - longitude
		@param cosl = cos(lambda)
		@param M - max order and harmonic
		@param g - output vector
	*/
	void earth_gravity(double r, double sint, double cost, double sinl, double cosl, int M, double* g)
	{
		double R = EARTH_RADIUS / r;
		double mu = GM / r;

		double xr1, xr2, xt1, xt2, xl1, xl2;
		double cnm, snm;
		double Pt;
		double* omega_r = new double[M + 1];
		double* omega_t = new double[M + 1];
		double* omega_l = new double[M + 1];

		double vr = 0.0, vt = 0.0, vl = 0.0;

		for (int m = 0; m < M + 1; m++)
		{
			xr1 = 0.0; xr2 = 0.0; xt1 = 0.0; xt2 = 0.0; xl1 = 0.0; xl2 = 0.0;

			for (int n = std::max(2, m); n < M + 1; n++)
			{
				if ((n == m) or (n - 1 < m)) Pt = (double)m * cost / sint * Ps[n][m];
				else Pt = ((double)n * cost * Ps[n][m] - f(n, m) * Ps[n - 1][m]) / sint;

				int ind = Order2[n][m];

				cnm = Cnm[ind] * pow(R, n);
				snm = Snm[ind] * pow(R, n);

				xl1 += cnm * Ps[n][m];
				xl2 += snm * Ps[n][m];

				xr1 -= (double)(n + 1) / r * cnm * Ps[n][m];
				xr2 -= (double)(n + 1) / r * snm * Ps[n][m];

				xt1 += cnm * Pt;
				xt2 += snm * Pt;
			}

			omega_r[m] = mu * (ml[1][m] * xr1 + ml[0][m] * xr2);
			omega_t[m] = mu * (ml[1][m] * xt1 + ml[0][m] * xt2);
			omega_l[m] = mu * (double)m * (ml[1][m] * xl2 - ml[0][m] * xl1);
		}

		// Modified forward column method

		double u;
		for (int m = M; m > 0; m--)
		{
			u = sqrt((double)(2 * m + 1) / (double)(2 * m)) * sint;
			if (m == 1) u = sqrt(3.0) * sint;

			vr = (vr + omega_r[m]) * u;
			vt = (vt + omega_t[m]) * u;
			vl = (vl + omega_l[m]) * u;
		}

		vr += omega_r[0];
		vt += omega_t[0];
		vl += omega_l[0];

		// for the central field
		vr -= mu / r;

		// components of gravitational acceleration
		g[0] = sint * cosl * vr + cost * cosl / r * vt - sinl * vl / r / sint;
		g[1] = sint * sinl * vr + cost * sinl / r * vt + cosl * vl / r / sint;
		g[2] = cost * vr - sint * vt / r;

		delete[] omega_r;
		delete[] omega_t;
		delete[] omega_l;
	}
}