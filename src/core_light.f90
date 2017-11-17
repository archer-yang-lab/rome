! --------------------------------------------------
SUBROUTINE core_light (delta, Kmat, Umat, Dvec, &
   & nobs, y, nlam, ulam, eps, &
    & maxit, anlam, npass, jerr, alpmat)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nlam
      INTEGER :: anlam
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: npass (nlam)
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: Kmat (nobs, nobs)
      DOUBLE PRECISION :: Umat (nobs, nobs)     
      DOUBLE PRECISION :: Dvec (nobs)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: alpmat (nobs+1, nlam)
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
     DOUBLE PRECISION :: r (nobs)
     DOUBLE PRECISION :: phi (nobs)
     DOUBLE PRECISION :: alpvec (nobs+1)
     DOUBLE PRECISION :: oalpvec (nobs+1)
     DOUBLE PRECISION :: dif (nobs+1)
     DOUBLE PRECISION :: Omega_lam_inv (nobs)
     DOUBLE PRECISION :: PIvec (nobs)
     DOUBLE PRECISION :: ONEvec (nobs)
     DOUBLE PRECISION :: Gscalar
     DOUBLE PRECISION :: tmp_vec1 (nobs)
     DOUBLE PRECISION :: tmp_vec2 (nobs)
     DOUBLE PRECISION :: tmp_scalar2
     DOUBLE PRECISION :: comb_tmp_vec1 (nobs+1)
     DOUBLE PRECISION :: comb_tmp_vec2 (nobs+1)
     DOUBLE PRECISION :: Ki (nobs+1, nobs)

    ! - - - begin - - -
   ONEvec = 1.0D0
   npass = 0
   r = y
   alpmat = 0.0D0
   alpvec = 0.0D0
   lambda_loop: DO l = 1,nlam
   dif = 0.0D0
   Ki = 1.0D0
   Ki(2:(nobs+1),:) = Kmat
   ! - - - computing Ku inverse - - - 
      Omega_lam_inv = 1.0D0/(Dvec*Dvec + ulam(l)*Dvec)
      PIvec = MATMUL(Umat,Omega_lam_inv*Dvec*MATMUL(TRANSPOSE(Umat), ONEvec))
      Gscalar = 1/(nobs-SUM(PIvec))
      ! update alpha
      alpha_loop: DO
         DO j = 1, nobs
                IF (ABS(r(j)) <= delta) THEN
                    phi (j) = r(j)
                ELSE
                    phi (j) = delta * SIGN(1.0D0, r(j))
                END IF
         END DO
         tmp_vec1 = phi-2.0D0*ulam(l)*alpvec(2:(nobs+1))
         tmp_scalar2 = SUM(phi) - DOT_PRODUCT(PIvec,MATMUL(Kmat, tmp_vec1))
         tmp_vec2 = MATMUL(Umat,Omega_lam_inv*Dvec*MATMUL(TRANSPOSE(Umat), tmp_vec1))
         comb_tmp_vec1(1) = 1.0D0
         comb_tmp_vec1(2:(nobs+1)) = -PIvec
         comb_tmp_vec2(1) = 0.0D0
         comb_tmp_vec2(2:(nobs+1)) = tmp_vec2
         oalpvec = alpvec
         alpvec =  oalpvec + 0.5D0*(tmp_scalar2*Gscalar*comb_tmp_vec1+comb_tmp_vec2)
         dif = alpvec - oalpvec
         r = r - MATMUL(dif, Ki)
         IF (DOT_PRODUCT(dif,dif)/DOT_PRODUCT(oalpvec,oalpvec) < eps) EXIT
         npass(l) = npass(l) + 1
         IF (SUM(npass) > maxit) EXIT
      ENDDO alpha_loop
      alpmat(:, l) = alpvec
      IF (SUM(npass) > maxit) THEN
         jerr = -l
         EXIT
      ENDIF
      anlam = l
   ENDDO lambda_loop

END SUBROUTINE core_light