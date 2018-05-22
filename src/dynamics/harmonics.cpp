//------------------------------------------------------------------------------
void Harmonic::CalculateField1(const Real& jday,  const Real pos[3], const Integer& nn,
                              const Integer& mm, const bool& fillgradient,
                              Real  acc[3],      Rmatrix33& gradient) const
{

   // calculate vector components ----------------------------------
   Real r = sqrt (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);    // Naming scheme from ref [3]
   Real s = pos[0]/r;
   Real t = pos[1]/r;
   Real u = pos[2]/r; // sin(phi), phi = geocentric latitude

   for (Integer n=0;  n<=NN;  ++n)
    for (Integer m=0;  m<=n && m<=MM;  ++m){
	   VR01[n][m] = sqrt(Real((n-m)*(n+m+1)));
     VR11[n][m] = sqrt(Real((2*n+1)*(n+m+2)*(n+m+1))/Real((2*n+3)));
	   if (m==0){
		   VR01[n][m] /= sqrt(Real(2));
		   VR11[n][m] /= sqrt(Real(2));
	   }
    }


   for (Integer m=0;  m<=MM+2;  ++m){
      for (Integer n=m+2;  n<=NN+2;  ++n)
      {
         N1[n][m] = sqrt (Real((2*n+1)*(2*n-1)) / Real((n-m)*(n+m)));
         N2[n][m] = sqrt (Real((2*n+1)*(n-m-1)*(n+m-1)) /
                          Real((2*n-3)*(n+m)*(n-m)));
      }
   }

   // initialize the diagonal elements (not a function of the input)
   A[0][0] = 1.0;
   for (Integer n=1;  n<=NN+2;  ++n)
      A[n][n] = sqrt (Real(2*n+1)/Real(2*n)) * A[n-1][n-1];

   // Calculate values for A -----------------------------------------
   // generate the off-diagonal elements
   A[1][0] = u*sqrt(Real(3.0));
   for (Integer n=1;  n<=NN+1 && n<=nn+1;  ++n)
      A[n+1][n] = u*sqrt(Real(2*n+3))*A[n][n];

   // apply column-fill recursion formula (Table 2, Row I, Ref.[1])
   for (Integer m=0;  m<=MM+1 && m<=mm+1;  ++m){
      for (Integer n=m+2;  n<=NN+1 && n<=nn+1;  ++n){
         A[n][m] = u * N1[n][m] * A[n-1][m] - N2[n][m] * A[n-2][m];
      }
      Re[m] = m==0 ? 1 : s*Re[m-1] - t*Im[m-1]; // real part of (s + i*t)^m
      Im[m] = m==0 ? 0 : s*Im[m-1] + t*Re[m-1]; // imaginary part of (s + i*t)^m
   }

   Real rho = bodyRadius/r;
   Real rho_np1 = -factor/r * rho;   // rho(0) ,Ref[3], Eq 26 , factor = mu for gravity
   Real a1 = 0;
   Real a2 = 0;
   Real a3 = 0;
   Real a4 = 0;
   Real sqrt2 = sqrt (Real(2));
   for (Integer n=1;  n<=NN && n<=nn;  ++n) {
      rho_np1 *= rho;
      Real sum1 = 0;
      Real sum2 = 0;
      Real sum3 = 0;
      Real sum4 = 0;

      for (Integer m=0;  m <= n && m<=MM && m<=mm;  ++m) {
         Real Cval,Sval;
         Cval = Cnm (jday,n,m);
         Sval = Snm (jday,n,m);
         // Pines Equation 27 (Part of)
         Real D =            (Cval*Re[m]   + Sval*Im[m]) * sqrt2;
         Real E = m==0 ? 0 : (Cval*Re[m-1] + Sval*Im[m-1]) * sqrt2;
         Real F = m==0 ? 0 : (Sval*Re[m-1] - Cval*Im[m-1]) * sqrt2;
         // Correct for normalization
         sum2 += m * A[n][m] * F;
         sum1 += m * A[n][m] * E;
         sum3 += VR01[n][m] * A[n][m+1] * D;
         sum4 += VR11[n][m] * A[n+1][m+1] * D;
      }
      // Pines Equation 30 and 30b (Part of)
      Real rr = rho_np1/bodyRadius;
      a1 += rr*sum1;
      a2 += rr*sum2;
      a3 += rr*sum3;
      a4 -= rr*sum4;
   }

   // Pines Equation 31
   acc[0] = a1+a4*s;
   acc[1] = a2+a4*t;
   acc[2] = a3+a4*u;

}
