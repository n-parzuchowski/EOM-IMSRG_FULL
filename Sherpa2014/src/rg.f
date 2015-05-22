      subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)

      integer n,nm,is1,is2,ierr,matz
      double precision a(nm,nm),wr(nm),wi(nm),z(nm,nm),fv1(nm)
      integer iv1(n)

c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real general matrix.

c     on input

c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.

c        n  is the order of the matrix  a.

c        a  contains the real general matrix.

c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.

c     on output

c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.  complex conjugate
c        pairs of eigenvalues appear consecutively with the
c        eigenvalue having the positive imaginary part first.

c        z  contains the real and imaginary parts of the eigenvectors
c        if matz is not zero.  if the j-th eigenvalue is real, the
c        j-th column of  z  contains its eigenvector.  if the j-th
c        eigenvalue is complex with positive imaginary part, the
c        j-th and (j+1)-th columns of  z  contain the real and
c        imaginary parts of its eigenvector.  the conjugate of this
c        vector is the eigenvector for the conjugate eigenvalue.

c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for hqr
c           and hqr2.  the normal completion code is zero.

c        iv1  and  fv1  are temporary storage arrays.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50

   10 call  balanc(nm,n,a,is1,is2,fv1)
      call  elmhes(nm,n,is1,is2,a,iv1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if (ierr .ne. 0) go to 50
      call  balbak(nm,n,is1,is2,fv1,n,z)
   50 return
      end
      subroutine balbak(nm,n,low,igh,scale,m,z)

      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),z(nm,m)
      double precision s

c     this subroutine is a translation of the algol procedure balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  balanc.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by  balanc.

c        scale contains information determining the permutations
c          and scaling factors used by  balanc.

c        m is the number of columns of z to be back transformed.

c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.

c     on output

c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120

      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s

  110 continue
c     ......... for i=low-1 step -1 until 1,
c               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140

         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue

  140 continue

  200 return
      end
      subroutine elmhes(nm,n,low,igh,a,int)

      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      double precision a(nm,n)
      double precision x,y
      integer int(igh)

c     this subroutine is a translation of the algol procedure elmhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.

c        a contains the input matrix.

c     on output

c        a contains the hessenberg matrix.  the multipliers
c          which were used in the reduction are stored in the
c          remaining triangle under the hessenberg matrix.

c        int contains information on the rows and columns
c          interchanged in the reduction.
c          only elements low through igh are used.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200

      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0d0
         i = m

         do 100 j = m, igh
            if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue

         int(m) = i
         if (i .eq. m) go to 130
c     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue

         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c     .......... end interchange ..........
  130    if (x .eq. 0.0d0) go to 180
         mp1 = m + 1

         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0d0) go to 160
            y = y / x
            a(i,mm1) = y

            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)

            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)

  160    continue

  180 continue

  200 return
      end
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)

      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas

c     this subroutine is a translation of the algol procedure hqr,
c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).

c     this subroutine finds the eigenvalues of a real
c     upper hessenberg matrix by the qr method.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.

c        h contains the upper hessenberg matrix.  information about
c          the transformations used in the reduction to hessenberg
c          form by  elmhes  or  orthes, if performed, is stored
c          in the remaining triangle under the hessenberg matrix.

c     on output

c        h has been destroyed.  therefore, it must be saved
c          before calling  hqr  if subsequent calculation and
c          back transformation of eigenvectors is to be performed.

c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.

c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated september 1989.

c     ------------------------------------------------------------------

      ierr = 0
      norm = 0.0d0
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n

         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))

         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue

      en = igh
      t = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x

      do 120 i = low, en
  120 h(i,i) = h(i,i) - x

      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue

  150 mp2 = m + 2

      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue

         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue

         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue

  260 continue

      go to 70
c     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.0d0
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      x = x + t
      if (q .lt. 0.0d0) go to 320
c     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine eltran(nm,n,low,igh,a,int,z)

      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),z(nm,n)
      integer int(igh)

c     this subroutine is a translation of the algol procedure elmtrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

c     this subroutine accumulates the stabilized elementary
c     similarity transformations used in the reduction of a
c     real general matrix to upper hessenberg form by  elmhes.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.

c        a contains the multipliers which were used in the
c          reduction by  elmhes  in its lower triangle
c          below the subdiagonal.

c        int contains information on the rows and columns
c          interchanged in the reduction by  elmhes.
c          only elements low through igh are used.

c     on output

c        z contains the transformation matrix produced in the
c          reduction by  elmhes.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

c     .......... initialize z to identity matrix ..........
      do 80 j = 1, n

         do 60 i = 1, n
   60    z(i,j) = 0.0d0

         z(j,j) = 1.0d0
   80 continue

      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1

         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)

         i = int(mp)
         if (i .eq. mp) go to 140

         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0d0
  130    continue

         z(i,mp) = 1.0d0
  140 continue

  200 return
      end
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)

      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,
     x        igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas

c     this subroutine is a translation of the algol procedure hqr2,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

c     this subroutine finds the eigenvalues and eigenvectors
c     of a real upper hessenberg matrix by the qr method.  the
c     eigenvectors of a real general matrix can also be found
c     if  elmhes  and  eltran  or  orthes  and  ortran  have
c     been used to reduce this general matrix to hessenberg form
c     and to accumulate the similarity transformations.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.

c        h contains the upper hessenberg matrix.

c        z contains the transformation matrix produced by  eltran
c          after the reduction by  elmhes, or by  ortran  after the
c          reduction by  orthes, if performed.  if the eigenvectors
c          of the hessenberg matrix are desired, z must contain the
c          identity matrix.

c     on output

c        h has been destroyed.

c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.

c        z contains the real and imaginary parts of the eigenvectors.
c          if the i-th eigenvalue is real, the i-th column of z
c          contains its eigenvector.  if the i-th eigenvalue is complex
c          with positive imaginary part, the i-th and (i+1)-th
c          columns of z contain the real and imaginary parts of its
c          eigenvector.  the eigenvectors are unnormalized.  if an
c          error exit is made, none of the eigenvectors has been found.

c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.

c     calls cdiv for complex division.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      ierr = 0
      norm = 0.0d0
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n

         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))

         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue

      en = igh
      t = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x

      do 120 i = low, en
  120 h(i,i) = h(i,i) - x

      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue

  150 mp2 = m + 2

      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue

         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
c     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue

         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
c     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue

  260 continue

      go to 70
c     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0d0
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0d0) go to 320
c     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = dabs(x) + dabs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(p*p+q*q)
      p = p / r
      q = q / r
c     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
c     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
c     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue

      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
c     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0d0
         if (na .eq. 0) go to 800
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0d0

            do 610 j = m, en
  610       r = r + h(i,j) * h(j,en)

            if (wi(i) .ge. 0.0d0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0d0) go to 640
            t = w
            if (t .ne. 0.0d0) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01d0 * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
c     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (dabs(x) .le. dabs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz

c     .......... overflow control ..........
  680       t = dabs(h(i,en))
            if (t .eq. 0.0d0) go to 700
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue

  700    continue
c     .......... end real vector ..........
         go to 800
c     .......... complex vector ..........
  710    m = na
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         if (dabs(h(en,na)) .le. dabs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cxdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0d0
            sa = 0.0d0

            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue

            if (wi(i) .ge. 0.0d0) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0d0) go to 780
            call cxdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
c     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0d0 * q
            if (vr .ne. 0.0d0 .or. vi .ne. 0.0d0) go to 784
               tst1 = norm * (dabs(w) + dabs(q) + dabs(x)
     x                      + dabs(y) + dabs(zz))
               vr = tst1
  783          vr = 0.01d0 * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cxdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,
     x                h(i,na),h(i,en))
            if (dabs(x) .le. dabs(zz) + dabs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cxdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,
     x                h(i+1,na),h(i+1,en))

c     .......... overflow control ..........
  790       t = dmax1(dabs(h(i,na)), dabs(h(i,en)))
            if (t .eq. 0.0d0) go to 795
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue

  795    continue
c     .......... end complex vector ..........
  800 continue
c     .......... end back substitution.
c                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840

         do 820 j = i, n
  820    z(i,j) = h(i,j)

  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)

         do 880 i = low, igh
            zz = 0.0d0

            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)

            z(i,j) = zz
  880 continue

      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine balanc(nm,n,a,low,igh,scale)

      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision a(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv

c     this subroutine is a translation of the algol procedure balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

c     this subroutine balances a real matrix and isolates
c     eigenvalues whenever possible.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        a contains the input matrix to be balanced.

c     on output

c        a contains the balanced matrix.

c        low and igh are two integers such that a(i,j)
c          is equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.

c        scale contains information determining the
c           permutations and scaling factors used.

c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j),      j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.

c     note that 1 is returned for igh if igh is zero formally.

c     the algol procedure exc contained in balance appears in
c     balanc  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      radix = 16.0d0

      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50

      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue

      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue

   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj

         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.0d0) go to 120
  110    continue

         m = l
         iexc = 1
         go to 20
  120 continue

      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1

  140 do 170 j = k, l

         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0d0) go to 170
  150    continue

         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.

      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0

         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(a(j,i))
            r = r + dabs(a(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.

         do 250 j = k, n
  250    a(i,j) = a(i,j) * g

         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f

  270 continue

      if (noconv) go to 190

  280 low = k
      igh = l
      return
      end

*     These routines are support routines only, provided since
*     complex arithmetic is not portable in FORTRAN77.  CXDIV and
*     CXMULT are complex division and multiplication; CXABS returns
*     the modulus of the supplied number; ICXMAX returns the index of
*     the complex entry of largest modulus in the arrays XRE and XIM;
*     CXCOPY returns a copy of the supplied complex number; CXAXPY
*     is a complex version of SAXPY; CXNRM2 returns the two-norm of
*     the supplied complex array; CXDOTU sums the term-by-term
*     product of the two supplied complex vectors; CXDOTC computes
*     the complex inner-product of the two supplied vectors; CXSCL
*     scales the supplied complex vector (XRE,XIM) by the complex
*     number (ARE,AIM); CXDSCL scales (XRE,XIM) by the double
*     precision number A; CXSQRT is the complex square-root;
*     CXASUM computes the sum of the moduli of the supplied complex
*     vector.

      SUBROUTINE CXDIV(AR,AI,BR,BI,CR,CI)
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI

*     complex division, (cr,ci) = (ar,ai)/(br,bi)

      DOUBLE PRECISION S,ARS,AIS,BRS,BIS
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END

      SUBROUTINE CXMULT(AR,AI,BR,BI,CR,CI)
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI

*     complex multiplication, (cr,ci) = (ar,ai)*(br,bi)

      CR = AR * BR - AI * BI
      CI = AR * BI + AI * BR
      RETURN
      END

      DOUBLE PRECISION FUNCTION CXABS(AR,AI)
      DOUBLE PRECISION AR, AI

*     complex absolute value

      DOUBLE PRECISION TRE,TIM,ANS,TMP
      TRE = ABS(AR)
      TIM = ABS(AI)
      IF (TRE .EQ. 0.0D0) THEN
         ANS = TIM
      ELSE IF (TIM .EQ. 0.0D0) THEN
         ANS = TRE
      ELSE IF (TRE .GT. TIM) THEN
         TMP = TIM/TRE
         ANS = TRE * SQRT(1.0D0 + TMP*TMP)
      ELSE
         TMP = TRE/TIM
         ANS = TIM * SQRT(1.0D0 + TMP*TMP)
      ENDIF
      CXABS = ANS
      RETURN
      END

      INTEGER FUNCTION ICXMAX(N,XRE,XIM,INCX)
      INTEGER N,INCX
      DOUBLE PRECISION XRE(*),XIM(*)

*     complex ICAMAX on (xre,xim)

      INTEGER I,IX
      DOUBLE PRECISION SMAX

      ICXMAX = 0
      IF( N.LT.1 ) RETURN
      ICXMAX = 1
      IF( N.EQ.1 ) RETURN
      IF(INCX.EQ.1) GO TO 30

*        code for increment not equal to 1

      IX = 1
      SMAX = ABS(XRE(1)) + ABS(XIM(1))
      IX = IX + INCX
      DO 20 I = 2,N
         IF(ABS(XRE(IX)) + ABS(XIM(IX)) .LE. SMAX) GO TO 10
         ICXMAX = I
         SMAX = ABS(XRE(IX)) + ABS(XIM(IX))
 10      IX = IX + INCX
 20   CONTINUE
      RETURN

*        code for increment equal to 1

 30   SMAX = ABS(XRE(1)) + ABS(XIM(1))
      DO 40 I = 2,N
         IF(ABS(XRE(I)) + ABS(XIM(I)) .LE. SMAX) GO TO 40
         ICXMAX = I
         SMAX = ABS(XRE(I))+ABS(XIM(I))
 40   CONTINUE
      RETURN
      END

      SUBROUTINE CXCOPY(N,XRE,XIM,YRE,YIM)
      INTEGER N
      DOUBLE PRECISION XRE(*),XIM(*),YRE(*),YIM(*)

*     complex copy on (yre,yim) <- (xre,xim)

      INTEGER I

      DO 10 I=1,N
         YRE(I) = XRE(I)
         YIM(I) = XIM(I)
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CXAXPY(N,ARE,AIM,XRE,XIM,YRE,YIM)
      INTEGER N
      DOUBLE PRECISION ARE,AIM,XRE(*),XIM(*),YRE(*),YIM(*)

*     complex axpy on (yre,yim) = (are,aim)*(xre,xim) + (yre,yim)

      INTEGER I

      DO 10 I=1,N
         YRE(I) = YRE(I) + ARE*XRE(I) - AIM*XIM(I)
         YIM(I) = YIM(I) + ARE*XIM(I) + AIM*XRE(I)
 10   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION CXNRM2(N,XRE,XIM)
      INTEGER N
      DOUBLE PRECISION XRE(*),XIM(*)

*     complex 2-norm on (xre,xim)

*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0)
      INTEGER I
      DOUBLE PRECISION SUM

      SUM = ZERO
      DO 10 I=1,N
         SUM = SUM + XRE(I)**2 + XIM(I)**2
 10   CONTINUE
      CXNRM2 = SQRT(SUM)
      RETURN
      END

      SUBROUTINE CXDOTU(N,XRE,XIM,YRE,YIM,ZRE,ZIM)
      INTEGER N
      DOUBLE PRECISION XRE(*),XIM(*),YRE(*),YIM(*),ZRE,ZIM

*     complex inner product (zre,zim) = (xre,xim)*(yre,yim)

*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0)
      INTEGER I

      ZRE = ZERO
      ZIM = ZERO
      DO 10 I=1,N
         ZRE = ZRE + XRE(I)*YRE(I) - XIM(I)*YIM(I)
         ZIM = ZIM + XRE(I)*YIM(I) + XIM(I)*YRE(I)
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CXDOTC(N,XRE,XIM,YRE,YIM,ZRE,ZIM)
      INTEGER N
      DOUBLE PRECISION XRE(*),XIM(*),YRE(*),YIM(*),ZRE,ZIM

*     complex conjugate inner product (zre,zim) = (xre,xim)*(yre,yim)

*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0)
      INTEGER I

      ZRE = ZERO
      ZIM = ZERO
      DO 10 I=1,N
         ZRE = ZRE + XRE(I)*YRE(I) + XIM(I)*YIM(I)
         ZIM = ZIM + XRE(I)*YIM(I) - XIM(I)*YRE(I)
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CXSCAL(N,ARE,AIM,XRE,XIM)
      INTEGER N
      DOUBLE PRECISION ARE,AIM,XRE(*),XIM(*)

*     complex scale (xre,xim) = (are,aim)*(xre,xim)

      INTEGER I
      DOUBLE PRECISION TRE,TIM

      DO 10 I=1,N
         TRE = ARE*XRE(I) - AIM*XIM(I)
         TIM = ARE*XIM(I) + AIM*XRE(I)
         XRE(I) = TRE
         XIM(I) = TIM
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CXDSCAL(N,A,XRE,XIM)
      INTEGER N
      DOUBLE PRECISION A,XRE(*),XIM(*)

*     complex scale with real (xre,xim) = a*(xre,xim)

      INTEGER I

      DO 10 I=1,N
         XRE(I) = A * XRE(I)
         XIM(I) = A * XIM(I)
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CXSQRT(XRE, XIM, YRE, YIM)
      DOUBLE PRECISION XRE,XIM,YRE,YIM

*     complex sqrt

      DOUBLE PRECISION X,Y,W,R
*     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,HALF
      PARAMETER (ZERO = 0.0,ONE = 1.0, TWO = 2.0 ,HALF = 0.5)
      IF ((XRE .EQ. ZERO) .AND. (XIM .EQ. ZERO)) THEN
         YRE = ZERO
         YIM = ZERO
         RETURN
      ELSE
         X = ABS(XRE)
         Y = ABS(XIM)
         IF (X .GE. Y) THEN
            R = Y / X
            W = SQRT(X)*SQRT(HALF*(ONE+SQRT(ONE+R*R)))
         ELSE
            R = X / Y
            W = SQRT(Y)*SQRT(HALF*(R+SQRT(ONE+R*R)))
         ENDIF
         IF (XRE .GE. ZERO) THEN
            YRE = W
            YIM = XIM / (TWO*W)
         ELSE
            IF (XIM .GE. ZERO) THEN
               YIM = W
            ELSE
               YIM = -W
            ENDIF
            YRE = XIM / (TWO* YIM)
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DCXASUM(N,XRE,XIM)
      INTEGER N
      DOUBLE PRECISION XRE(*),XIM(*)
*
*     compute complex sum of absolute values in (xre,xim)

*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0)
      INTEGER I
      DOUBLE PRECISION SUM

      SUM = ZERO
      DO 10 I=1,N
         SUM = SUM + ABS(XRE(I)) + ABS(XIM(I))
 10   CONTINUE
      DCXASUM = SUM
      RETURN
      END

      SUBROUTINE CXSIGN1(XRE,XIM,YRE,YIM,ZRE,ZIM)
      DOUBLE PRECISION XRE,XIM,YRE,YIM,ZRE,ZIM

      DOUBLE PRECISION TMP1, TMP2

      TMP1 = ABS(XRE) + ABS(XIM)
      TMP2 = ABS(YRE) + ABS(YIM)
      ZRE = (YRE / TMP2) * TMP1
      ZIM = (YIM / TMP2) * TMP1
      RETURN
      END
