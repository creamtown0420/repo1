subroutine dqag_sk3d(f, xa, xb, ya, yb, za, zb, eps, s, ier, key)
   implicit none
   integer, intent(in)::key
   real*16, intent(in)::xa, xb, ya, yb, za, zb, eps
   real*16, intent(out)::s
   integer, intent(out)::ier
   ! ier = 0 : success, converged.
   ! ier > 0 : fail, not converged.

   integer::neval, limit, lenw, last
   integer, allocatable::iwork(:)
   real*16, allocatable::work(:)
   real*16::epsabs, abserr
   interface
      function f(x, y, z)
         implicit none
         real*16, intent(in)::x, y, z
         real*16::f
      end function f
   end interface
   ier = 1
   epsabs = -1q0
   !limit=200
   limit = 1000
   lenw = limit*4
   allocate (iwork(1:limit), work(1:lenw))

   s = 0q0
   call dqag3(f, xa, xb, epsabs, eps, key, s, abserr, neval, ier, &
              limit, lenw, last, iwork, work, ya, yb, za, zb)

   deallocate (iwork, work)
   return
end subroutine dqag_sk3d

subroutine dqag3(f, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
                 limit, lenw, last, iwork, work, ya, yb, za, zb)
   implicit none
   real*16 a, abserr, b, epsabs, epsrel, result, work, ya, yb, za, zb
   integer ier, iwork, key, last, lenw, limit, lvl, l1, l2, l3, neval
   dimension iwork(limit), work(lenw)
   interface
      function f(x, y, z)
         implicit none
         real*16, intent(in)::x, y, z
         real*16::f
      end function f
   end interface
   ier = 6
   neval = 0
   last = 0
   result = 0.0q0 !0.0q0
   abserr = 0.0q0
   if (limit .ge. 1 .and. lenw .ge. limit*4) then
      l1 = limit + 1
      l2 = limit + l1
      l3 = limit + l2
      call dqage3(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
                  ier, work(1), work(l1), work(l2), work(l3), iwork, last, ya, yb, za, zb)
      lvl = 0
   end if
   if (ier .eq. 6) lvl = 1
   if (ier .ne. 0) call xerror("26habnormal return from dqag", 26, ier, lvl)

   return
end subroutine dqag3

subroutine dqage3(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                  neval, ier, alist, blist, rlist, elist, iord, last, ya, yb, za, zb)
   implicit none
   real*16 a, abserr, alist, area, area1, area12, area2, a1, a2, b, &
      blist, b1, b2, abs, defabs, defab1, defab2, dmax1, elist, epmach, &
      epsabs, epsrel, errbnd, errmax, error1, error2, erro12, errsum, &
      resabs, result, rlist, uflow, ya, yb, za, zb
   integer ier, iord, iroff1, iroff2, k, key, keyf, last, limit, maxerr, neval, &
      nrmax, igt, igk
   dimension alist(limit), blist(limit), elist(limit), iord(limit), &
      rlist(limit)
   interface
      function f(x, y, z)
         implicit none
         real*16, intent(in)::x, y, z
         real*16::f
      end function f
   end interface
   epmach = epsilon(a)
   uflow = tiny(a)
   ier = 0
   neval = 0
   last = 0
   result = 0.0q0
   abserr = 0.0q0
   alist(1) = a
   blist(1) = b
   rlist(1) = 0.0q0
   elist(1) = 0.0q0
   iord(1) = 0
   if (epsabs .le. 0.0q0 .and. &
       epsrel .lt. dmax1(0.5d+02*epmach, 0.5d-28)) ier = 6
   if (ier .ne. 6) then
      keyf = key
      if (key .le. 0) keyf = 1
      if (key .ge. 7) keyf = 6
      neval = 0
      if (keyf .eq. 1) call dqk15_3(f, a, b, result, abserr, defabs, resabs, ya, yb, za, zb, epsrel)
      if (keyf .eq. 2) call dqk21_3(f, a, b, result, abserr, defabs, resabs, ya, yb, za, zb, epsrel)
      if (keyf .eq. 3) call dqk31_3(f, a, b, result, abserr, defabs, resabs, ya, yb, za, zb, epsrel)
      if (keyf .eq. 4) call dqk41_3(f, a, b, result, abserr, defabs, resabs, ya, yb, za, zb, epsrel)
      if (keyf .eq. 5) call dqk51_3(f, a, b, result, abserr, defabs, resabs, ya, yb, za, zb, epsrel)
      if (keyf .eq. 6) call dqk61_3(f, a, b, result, abserr, defabs, resabs, ya, yb, za, zb, epsrel)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      errbnd = dmax1(epsabs, epsrel*abs(result))
      if (abserr .le. 0.5d+02*epmach*defabs .and. abserr .gt. errbnd) ier = 2
      if (limit .eq. 1) ier = 1
      igk = 0
      if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr .ne. resabs) &
          .or. abserr .eq. 0.0q0) igk = 60
      if (igk .ne. 60) then
         errmax = abserr
         maxerr = 1
         area = result
         errsum = abserr
         nrmax = 1
         iroff1 = 0
         iroff2 = 0
         do last = 2, limit
            a1 = alist(maxerr)
            b1 = 0.5q0*(alist(maxerr) + blist(maxerr))
            a2 = b1
            b2 = blist(maxerr)
            if (keyf .eq. 1) call dqk15_3(f, a1, b1, area1, error1, resabs, defab1, ya, yb, za, zb, epsrel)
            if (keyf .eq. 2) call dqk21_3(f, a1, b1, area1, error1, resabs, defab1, ya, yb, za, zb, epsrel)
            if (keyf .eq. 3) call dqk31_3(f, a1, b1, area1, error1, resabs, defab1, ya, yb, za, zb, epsrel)
            if (keyf .eq. 4) call dqk41_3(f, a1, b1, area1, error1, resabs, defab1, ya, yb, za, zb, epsrel)
            if (keyf .eq. 5) call dqk51_3(f, a1, b1, area1, error1, resabs, defab1, ya, yb, za, zb, epsrel)
            if (keyf .eq. 6) call dqk61_3(f, a1, b1, area1, error1, resabs, defab1, ya, yb, za, zb, epsrel)
            if (keyf .eq. 1) call dqk15_3(f, a2, b2, area2, error2, resabs, defab2, ya, yb, za, zb, epsrel)
            if (keyf .eq. 2) call dqk21_3(f, a2, b2, area2, error2, resabs, defab2, ya, yb, za, zb, epsrel)
            if (keyf .eq. 3) call dqk31_3(f, a2, b2, area2, error2, resabs, defab2, ya, yb, za, zb, epsrel)
            if (keyf .eq. 4) call dqk41_3(f, a2, b2, area2, error2, resabs, defab2, ya, yb, za, zb, epsrel)
            if (keyf .eq. 5) call dqk51_3(f, a2, b2, area2, error2, resabs, defab2, ya, yb, za, zb, epsrel)
            if (keyf .eq. 6) call dqk61_3(f, a2, b2, area2, error2, resabs, defab2, ya, yb, za, zb, epsrel)
            neval = neval + 1
            area12 = area1 + area2
            erro12 = error1 + error2
            errsum = errsum + erro12 - errmax
            area = area + area12 - rlist(maxerr)
            if (defab1 .ne. error1 .and. defab2 .ne. error2) then
               if (abs(rlist(maxerr) - area12) .le. 0.1d-04*abs(area12) &
                   .and. erro12 .ge. 0.99q0*errmax) iroff1 = iroff1 + 1
               if (last .gt. 10 .and. erro12 .gt. errmax) iroff2 = iroff2 + 1
            end if
            rlist(maxerr) = area1
            rlist(last) = area2
            errbnd = dmax1(epsabs, epsrel*abs(area))
            if (errsum .gt. errbnd) then
               if (iroff1 .ge. 6 .or. iroff2 .ge. 20) ier = 2
               if (last .eq. limit) ier = 1
               if (dmax1(abs(a1), abs(b2)) .le. (0.1d+01 + 0.1d+03* &
                                                 epmach)*(abs(a2) + 0.1d+04*uflow)) ier = 3
            end if
            igt = 0
            if (error2 .le. error1) then
               alist(last) = a2
               blist(maxerr) = b1
               blist(last) = b2
               elist(maxerr) = error1
               elist(last) = error2
               igt = 20
            end if
            if (igt .ne. 20) then
               alist(maxerr) = a2
               alist(last) = a1
               blist(last) = b1
               rlist(maxerr) = area2
               rlist(last) = area1
               elist(maxerr) = error2
               elist(last) = error1
            end if
            call dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax)
            if (ier .ne. 0 .or. errsum .le. errbnd) exit
         end do
         if (last .eq. limit + 1) last = last - 1
         result = 0.0q0
         do k = 1, last
            result = result + rlist(k)
         end do
         abserr = errsum
      end if
      if (keyf .ne. 1) neval = (10*keyf + 1)*(2*neval + 1)
      if (keyf .eq. 1) neval = 30*neval + 15
   end if
   return
end subroutine dqage3

subroutine dqk15_3(f, a, b, result, abserr, resabs, resasc, ya, yb, za, zb, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ya, yb, za, zb, eps
   integer j, jtw, jtwm1, ier
   interface
      function f(x, y, z)
         implicit none
         real*16, intent(in)::x, y, z
         real*16::f
      end function f
   end interface
   dimension fv1(7), fv2(7), wg(4), wgk(8), xgk(8)
   data wg(1)/0.129484966168869693270611432679082q0/
   data wg(2)/0.279705391489276667901467771423780q0/
   data wg(3)/0.381830050505118944950369775488975q0/
   data wg(4)/0.417959183673469387755102040816327q0/

   data xgk(1)/0.991455371120812639206854697526329q0/
   data xgk(2)/0.949107912342758524526189684047851q0/
   data xgk(3)/0.864864423359769072789712788640926q0/
   data xgk(4)/0.741531185599394439863864773280788q0/
   data xgk(5)/0.586087235467691130294144838258730q0/
   data xgk(6)/0.405845151377397166906606412076961q0/
   data xgk(7)/0.207784955007898467600689403773245q0/
   data xgk(8)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.022935322010529224963732008058970q0/
   data wgk(2)/0.063092092629978553290700663189204q0/
   data wgk(3)/0.104790010322250183839876322541518q0/
   data wgk(4)/0.140653259715525918745189590510238q0/
   data wgk(5)/0.169004726639267902826583426598550q0/
   data wgk(6)/0.190350578064785409913256402421014q0/
   data wgk(7)/0.204432940075298892414161999234649q0/
   data wgk(8)/0.209482141084727828012999174891714q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   call dqag_sk3e(f, ya, yb, eps, fc, ier, 1, za, zb, centr)
   resg = fc*wg(4)
   resk = fc*wgk(8)
   resabs = abs(resk)
   do j = 1, 3
      jtw = j*2
      absc = hlgth*xgk(jtw)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 1, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 1, za, zb, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 4
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 1, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 1, za, zb, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(8)*abs(fc - reskh)
   do j = 1, 7
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk15_3
!===================================ここまでinterfaceやった
subroutine dqk21_3(f, a, b, result, abserr, resabs, resasc, ya, yb, za, zb, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ya, yb, za, zb, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(10), fv2(10), wg(5), wgk(11), xgk(11)

   data wg(1)/0.066671344308688137593568809893332q0/
   data wg(2)/0.149451349150580593145776339657697q0/
   data wg(3)/0.219086362515982043995534934228163q0/
   data wg(4)/0.269266719309996355091226921569469q0/
   data wg(5)/0.295524224714752870173892994651338q0/

   data xgk(1)/0.995657163025808080735527280689003q0/
   data xgk(2)/0.973906528517171720077964012084452q0/
   data xgk(3)/0.930157491355708226001207180059508q0/
   data xgk(4)/0.865063366688984510732096688423493q0/
   data xgk(5)/0.780817726586416897063717578345042q0/
   data xgk(6)/0.679409568299024406234327365114874q0/
   data xgk(7)/0.562757134668604683339000099272694q0/
   data xgk(8)/0.433395394129247190799265943165784q0/
   data xgk(9)/0.294392862701460198131126603103866q0/
   data xgk(10)/0.148874338981631210884826001129720q0/
   data xgk(11)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.011694638867371874278064396062192q0/
   data wgk(2)/0.032558162307964727478818972459390q0/
   data wgk(3)/0.054755896574351996031381300244580q0/
   data wgk(4)/0.075039674810919952767043140916190q0/
   data wgk(5)/0.093125454583697605535065465083366q0/
   data wgk(6)/0.109387158802297641899210590325805q0/
   data wgk(7)/0.123491976262065851077958109831074q0/
   data wgk(8)/0.134709217311473325928054001771707q0/
   data wgk(9)/0.142775938577060080797094273138717q0/
   data wgk(10)/0.147739104901338491374841515972068q0/
   data wgk(11)/0.149445554002916905664936468389821q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   call dqag_sk3e(f, ya, yb, eps, fc, ier, 2, za, zb, centr)
   resk = wgk(11)*fc
   resabs = abs(resk)
   do j = 1, 5
      jtw = 2*j
      absc = hlgth*xgk(jtw)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,2,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,2,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 2, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 2, za, zb, centr + absc)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 5
      jtwm1 = 2*j - 1
      absc = hlgth*xgk(jtwm1)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,2,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,2,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 2, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 2, za, zb, centr + absc)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(11)*abs(fc - reskh)
   do j = 1, 10
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk21_3

subroutine dqk31_3(f, a, b, result, abserr, resabs, resasc, ya, yb, za, zb, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ya, yb, za, zb, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(15), fv2(15), xgk(16), wgk(16), wg(8)

   data wg(1)/0.030753241996117268354628393577204q0/
   data wg(2)/0.070366047488108124709267416450667q0/
   data wg(3)/0.107159220467171935011869546685869q0/
   data wg(4)/0.139570677926154314447804794511028q0/
   data wg(5)/0.166269205816993933553200860481209q0/
   data wg(6)/0.186161000015562211026800561866423q0/
   data wg(7)/0.198431485327111576456118326443839q0/
   data wg(8)/0.202578241925561272880620199967519q0/

   data xgk(1)/0.998002298693397060285172840152271q0/
   data xgk(2)/0.987992518020485428489565718586613q0/
   data xgk(3)/0.967739075679139134257347978784337q0/
   data xgk(4)/0.937273392400705904307758947710209q0/
   data xgk(5)/0.897264532344081900882509656454496q0/
   data xgk(6)/0.848206583410427216200648320774217q0/
   data xgk(7)/0.790418501442465932967649294817947q0/
   data xgk(8)/0.724417731360170047416186054613938q0/
   data xgk(9)/0.650996741297416970533735895313275q0/
   data xgk(10)/0.570972172608538847537226737253911q0/
   data xgk(11)/0.485081863640239680693655740232351q0/
   data xgk(12)/0.394151347077563369897207370981045q0/
   data xgk(13)/0.299180007153168812166780024266389q0/
   data xgk(14)/0.201194093997434522300628303394596q0/
   data xgk(15)/0.101142066918717499027074231447392q0/
   data xgk(16)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.005377479872923348987792051430128q0/
   data wgk(2)/0.015007947329316122538374763075807q0/
   data wgk(3)/0.025460847326715320186874001019653q0/
   data wgk(4)/0.035346360791375846222037948478360q0/
   data wgk(5)/0.044589751324764876608227299373280q0/
   data wgk(6)/0.053481524690928087265343147239430q0/
   data wgk(7)/0.062009567800670640285139230960803q0/
   data wgk(8)/0.069854121318728258709520077099147q0/
   data wgk(9)/0.076849680757720378894432777482659q0/
   data wgk(10)/0.083080502823133021038289247286104q0/
   data wgk(11)/0.088564443056211770647275443693774q0/
   data wgk(12)/0.093126598170825321225486872747346q0/
   data wgk(13)/0.096642726983623678505179907627589q0/
   data wgk(14)/0.099173598721791959332393173484603q0/
   data wgk(15)/0.100769845523875595044946662617570q0/
   data wgk(16)/0.101330007014791549017374792767493q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   !call dqag_sk2e(f,ya,yb,eps,fc,ier,3,centr)
   call dqag_sk3e(f, ya, yb, eps, fc, ier, 3, za, zb, centr)
   resg = wg(8)*fc
   resk = wgk(16)*fc
   resabs = abs(resk)
   do j = 1, 7
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,3,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,3,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 3, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 3, za, zb, centr + absc)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 8
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,3,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,3,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 3, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 3, za, zb, centr + absc)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(16)*abs(fc - reskh)
   do j = 1, 15
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk31_3

subroutine dqk41_3(f, a, b, result, abserr, resabs, resasc, ya, yb, za, zb, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ya, yb, za, zb, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(20), fv2(20), xgk(21), wgk(21), wg(10)

   data wg(1)/0.017614007139152118311861962351853q0/
   data wg(2)/0.040601429800386941331039952274932q0/
   data wg(3)/0.062672048334109063569506535187042q0/
   data wg(4)/0.083276741576704748724758143222046q0/
   data wg(5)/0.101930119817240435036750135480350q0/
   data wg(6)/0.118194531961518417312377377711382q0/
   data wg(7)/0.131688638449176626898494499748163q0/
   data wg(8)/0.142096109318382051329298325067165q0/
   data wg(9)/0.149172986472603746787828737001969q0/
   data wg(10)/0.152753387130725850698084331955098q0/
   !
   data xgk(1)/0.998859031588277663838315576545863q0/
   data xgk(2)/0.993128599185094924786122388471320q0/
   data xgk(3)/0.981507877450250259193342994720217q0/
   data xgk(4)/0.963971927277913791267666131197277q0/
   data xgk(5)/0.940822633831754753519982722212443q0/
   data xgk(6)/0.912234428251325905867752441203298q0/
   data xgk(7)/0.878276811252281976077442995113078q0/
   data xgk(8)/0.839116971822218823394529061701521q0/
   data xgk(9)/0.795041428837551198350638833272788q0/
   data xgk(10)/0.746331906460150792614305070355642q0/
   data xgk(11)/0.693237656334751384805490711845932q0/
   data xgk(12)/0.636053680726515025452836696226286q0/
   data xgk(13)/0.575140446819710315342946036586425q0/
   data xgk(14)/0.510867001950827098004364050955251q0/
   data xgk(15)/0.443593175238725103199992213492640q0/
   data xgk(16)/0.373706088715419560672548177024927q0/
   data xgk(17)/0.301627868114913004320555356858592q0/
   data xgk(18)/0.227785851141645078080496195368575q0/
   data xgk(19)/0.152605465240922675505220241022678q0/
   data xgk(20)/0.076526521133497333754640409398838q0/
   data xgk(21)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.003073583718520531501218293246031q0/
   data wgk(2)/0.008600269855642942198661787950102q0/
   data wgk(3)/0.014626169256971252983787960308868q0/
   data wgk(4)/0.020388373461266523598010231432755q0/
   data wgk(5)/0.025882133604951158834505067096153q0/
   data wgk(6)/0.031287306777032798958543119323801q0/
   data wgk(7)/0.036600169758200798030557240707211q0/
   data wgk(8)/0.041668873327973686263788305936895q0/
   data wgk(9)/0.046434821867497674720231880926108q0/
   data wgk(10)/0.050944573923728691932707670050345q0/
   data wgk(11)/0.055195105348285994744832372419777q0/
   data wgk(12)/0.059111400880639572374967220648594q0/
   data wgk(13)/0.062653237554781168025870122174255q0/
   data wgk(14)/0.065834597133618422111563556969398q0/
   data wgk(15)/0.068648672928521619345623411885368q0/
   data wgk(16)/0.071054423553444068305790361723210q0/
   data wgk(17)/0.073030690332786667495189417658913q0/
   data wgk(18)/0.074582875400499188986581418362488q0/
   data wgk(19)/0.075704497684556674659542775376617q0/
   data wgk(20)/0.076377867672080736705502835038061q0/
   data wgk(21)/0.076600711917999656445049901530102q0/

   epmach = epsilon(a)
   uflow = tiny(a)

   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   !call dqag_sk2e(f,ya,yb,eps,fc,ier,4,centr)
   call dqag_sk3e(f, ya, yb, eps, fc, ier, 4, za, zb, centr)

   resk = wgk(21)*fc
   resabs = abs(resk)
   do j = 1, 10
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,4,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,4,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 4, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 4, za, zb, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 10
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,4,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,4,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 4, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 4, za, zb, centr + absc)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(21)*abs(fc - reskh)
   do j = 1, 20
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk41_3

subroutine dqk51_3(f, a, b, result, abserr, resabs, resasc, ya, yb, za, zb, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ya, yb, za, zb, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(25), fv2(25), xgk(26), wgk(26), wg(13)

   data wg(1)/0.011393798501026287947902964113235q0/
   data wg(2)/0.026354986615032137261901815295299q0/
   data wg(3)/0.040939156701306312655623487711646q0/
   data wg(4)/0.054904695975835191925936891540473q0/
   data wg(5)/0.068038333812356917207187185656708q0/
   data wg(6)/0.080140700335001018013234959669111q0/
   data wg(7)/0.091028261982963649811497220702892q0/
   data wg(8)/0.100535949067050644202206890392686q0/
   data wg(9)/0.108519624474263653116093957050117q0/
   data wg(10)/0.114858259145711648339325545869556q0/
   data wg(11)/0.119455763535784772228178126512901q0/
   data wg(12)/0.122242442990310041688959518945852q0/
   data wg(13)/0.123176053726715451203902873079050q0/
   !
   data xgk(1)/0.999262104992609834193457486540341q0/
   data xgk(2)/0.995556969790498097908784946893902q0/
   data xgk(3)/0.988035794534077247637331014577406q0/
   data xgk(4)/0.976663921459517511498315386479594q0/
   data xgk(5)/0.961614986425842512418130033660167q0/
   data xgk(6)/0.942974571228974339414011169658471q0/
   data xgk(7)/0.920747115281701561746346084546331q0/
   data xgk(8)/0.894991997878275368851042006782805q0/
   data xgk(9)/0.865847065293275595448996969588340q0/
   data xgk(10)/0.833442628760834001421021108693570q0/
   data xgk(11)/0.797873797998500059410410904994307q0/
   data xgk(12)/0.759259263037357630577282865204361q0/
   data xgk(13)/0.717766406813084388186654079773298q0/
   data xgk(14)/0.673566368473468364485120633247622q0/
   data xgk(15)/0.626810099010317412788122681624518q0/
   data xgk(16)/0.577662930241222967723689841612654q0/
   data xgk(17)/0.526325284334719182599623778158010q0/
   data xgk(18)/0.473002731445714960522182115009192q0/
   data xgk(19)/0.417885382193037748851814394594572q0/
   data xgk(20)/0.361172305809387837735821730127641q0/
   data xgk(21)/0.303089538931107830167478909980339q0/
   data xgk(22)/0.243866883720988432045190362797452q0/
   data xgk(23)/0.183718939421048892015969888759528q0/
   data xgk(24)/0.122864692610710396387359818808037q0/
   data xgk(25)/0.061544483005685078886546392366797q0/
   data xgk(26)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.001987383892330315926507851882843q0/
   data wgk(2)/0.005561932135356713758040236901066q0/
   data wgk(3)/0.009473973386174151607207710523655q0/
   data wgk(4)/0.013236229195571674813656405846976q0/
   data wgk(5)/0.016847817709128298231516667536336q0/
   data wgk(6)/0.020435371145882835456568292235939q0/
   data wgk(7)/0.024009945606953216220092489164881q0/
   data wgk(8)/0.027475317587851737802948455517811q0/
   data wgk(9)/0.030792300167387488891109020215229q0/
   data wgk(10)/0.034002130274329337836748795229551q0/
   data wgk(11)/0.037116271483415543560330625367620q0/
   data wgk(12)/0.040083825504032382074839284467076q0/
   data wgk(13)/0.042872845020170049476895792439495q0/
   data wgk(14)/0.045502913049921788909870584752660q0/
   data wgk(15)/0.047982537138836713906392255756915q0/
   data wgk(16)/0.050277679080715671963325259433440q0/
   data wgk(17)/0.052362885806407475864366712137873q0/
   data wgk(18)/0.054251129888545490144543370459876q0/
   data wgk(19)/0.055950811220412317308240686382747q0/
   data wgk(20)/0.057437116361567832853582693939506q0/
   data wgk(21)/0.058689680022394207961974175856788q0/
   data wgk(22)/0.059720340324174059979099291932562q0/
   data wgk(23)/0.060539455376045862945360267517565q0/
   data wgk(24)/0.061128509717053048305859030416293q0/
   data wgk(25)/0.061471189871425316661544131965264q0/
   !       note: wgk (26) was calculated from the values of wgk(1..25)
   data wgk(26)/0.061580818067832935078759824240066q0/

   epmach = epsilon(a)
   uflow = tiny(a)

   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)
   !  call dqag_sk2e(f,ya,yb,eps,fc,ier,5,centr)
   call dqag_sk3e(f, ya, yb, eps, fc, ier, 5, za, zb, centr)
   resg = wg(13)*fc
   resk = wgk(26)*fc
   resabs = abs(resk)
   do j = 1, 12
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,5,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,5,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 5, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 5, za, zb, centr + absc)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 13
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,5,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,5,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 5, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 5, za, zb, centr + absc)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(26)*abs(fc - reskh)
   do j = 1, 25
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk51_3

subroutine dqk61_3(f, a, b, result, abserr, resabs, resasc, ya, yb, za, zb, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ya, yb, za, zb, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(30), fv2(30), xgk(31), wgk(31), wg(15)

   data wg(1)/0.007968192496166605615465883474674q0/
   data wg(2)/0.018466468311090959142302131912047q0/
   data wg(3)/0.028784707883323369349719179611292q0/
   data wg(4)/0.038799192569627049596801936446348q0/
   data wg(5)/0.048402672830594052902938140422808q0/
   data wg(6)/0.057493156217619066481721689402056q0/
   data wg(7)/0.065974229882180495128128515115962q0/
   data wg(8)/0.073755974737705206268243850022191q0/
   data wg(9)/0.080755895229420215354694938460530q0/
   data wg(10)/0.086899787201082979802387530715126q0/
   data wg(11)/0.092122522237786128717632707087619q0/
   data wg(12)/0.096368737174644259639468626351810q0/
   data wg(13)/0.099593420586795267062780282103569q0/
   data wg(14)/0.101762389748405504596428952168554q0/
   data wg(15)/0.102852652893558840341285636705415q0/
   !
   data xgk(1)/0.999484410050490637571325895705811q0/
   data xgk(2)/0.996893484074649540271630050918695q0/
   data xgk(3)/0.991630996870404594858628366109486q0/
   data xgk(4)/0.983668123279747209970032581605663q0/
   data xgk(5)/0.973116322501126268374693868423707q0/
   data xgk(6)/0.960021864968307512216871025581798q0/
   data xgk(7)/0.944374444748559979415831324037439q0/
   data xgk(8)/0.926200047429274325879324277080474q0/
   data xgk(9)/0.905573307699907798546522558925958q0/
   data xgk(10)/0.882560535792052681543116462530226q0/
   data xgk(11)/0.857205233546061098958658510658944q0/
   data xgk(12)/0.829565762382768397442898119732502q0/
   data xgk(13)/0.799727835821839083013668942322683q0/
   data xgk(14)/0.767777432104826194917977340974503q0/
   data xgk(15)/0.733790062453226804726171131369528q0/
   data xgk(16)/0.697850494793315796932292388026640q0/
   data xgk(17)/0.660061064126626961370053668149271q0/
   data xgk(18)/0.620526182989242861140477556431189q0/
   data xgk(19)/0.579345235826361691756024932172540q0/
   data xgk(20)/0.536624148142019899264169793311073q0/
   data xgk(21)/0.492480467861778574993693061207709q0/
   data xgk(22)/0.447033769538089176780609900322854q0/
   data xgk(23)/0.400401254830394392535476211542661q0/
   data xgk(24)/0.352704725530878113471037207089374q0/
   data xgk(25)/0.304073202273625077372677107199257q0/
   data xgk(26)/0.254636926167889846439805129817805q0/
   data xgk(27)/0.204525116682309891438957671002025q0/
   data xgk(28)/0.153869913608583546963794672743256q0/
   data xgk(29)/0.102806937966737030147096751318001q0/
   data xgk(30)/0.051471842555317695833025213166723q0/
   data xgk(31)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.001389013698677007624551591226760q0/
   data wgk(2)/0.003890461127099884051267201844516q0/
   data wgk(3)/0.006630703915931292173319826369750q0/
   data wgk(4)/0.009273279659517763428441146892024q0/
   data wgk(5)/0.011823015253496341742232898853251q0/
   data wgk(6)/0.014369729507045804812451432443580q0/
   data wgk(7)/0.016920889189053272627572289420322q0/
   data wgk(8)/0.019414141193942381173408951050128q0/
   data wgk(9)/0.021828035821609192297167485738339q0/
   data wgk(10)/0.024191162078080601365686370725232q0/
   data wgk(11)/0.026509954882333101610601709335075q0/
   data wgk(12)/0.028754048765041292843978785354334q0/
   data wgk(13)/0.030907257562387762472884252943092q0/
   data wgk(14)/0.032981447057483726031814191016854q0/
   data wgk(15)/0.034979338028060024137499670731468q0/
   data wgk(16)/0.036882364651821229223911065617136q0/
   data wgk(17)/0.038678945624727592950348651532281q0/
   data wgk(18)/0.040374538951535959111995279752468q0/
   data wgk(19)/0.041969810215164246147147541285970q0/
   data wgk(20)/0.043452539701356069316831728117073q0/
   data wgk(21)/0.044814800133162663192355551616723q0/
   data wgk(22)/0.046059238271006988116271735559374q0/
   data wgk(23)/0.047185546569299153945261478181099q0/
   data wgk(24)/0.048185861757087129140779492298305q0/
   data wgk(25)/0.049055434555029778887528165367238q0/
   data wgk(26)/0.049795683427074206357811569379942q0/
   data wgk(27)/0.050405921402782346840893085653585q0/
   data wgk(28)/0.050881795898749606492297473049805q0/
   data wgk(29)/0.051221547849258772170656282604944q0/
   data wgk(30)/0.051426128537459025933862879215781q0/
   data wgk(31)/0.051494729429451567558340433647099q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   !
   centr = 0.5q0*(b + a)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   !call dqag_sk2e(f,ya,yb,eps,fc,ier,6,centr)
   call dqag_sk3e(f, ya, yb, eps, fc, ier, 6, za, zb, centr)

   resk = wgk(31)*fc
   resabs = abs(resk)
   do j = 1, 15
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,6,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,6,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 6, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 6, za, zb, centr + absc)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 15
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !call dqag_sk2e(f,ya,yb,eps,fval1,ier,6,centr-absc)
      !call dqag_sk2e(f,ya,yb,eps,fval2,ier,6,centr+absc)
      call dqag_sk3e(f, ya, yb, eps, fval1, ier, 6, za, zb, centr - absc)
      call dqag_sk3e(f, ya, yb, eps, fval2, ier, 6, za, zb, centr + absc)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(31)*abs(fc - reskh)
   do j = 1, 30
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk61_3

!//////

subroutine dqag_sk3e(f, a, b, eps, s, ier, key, za, zb, tx)
   implicit none
   integer, intent(in)::key
   real*16, intent(in)::a, b, eps, tx, za, zb
   real*16, intent(out)::s
   integer, intent(out)::ier
   real*16::f
   ! ier = 0 : success, converged.
   ! ier > 0 : fail, not converged.

   integer::neval, limit, lenw, last
   integer, allocatable::iwork(:)
   real*16, allocatable::work(:)
   real*16::epsabs, abserr
   external::f

   ier = 1
   epsabs = -1q0
   !limit=200
   limit = 1000
   lenw = limit*4
   allocate (iwork(1:limit), work(1:lenw))

   s = 0q0
   call dqag3e(f, a, b, epsabs, eps, key, s, abserr, neval, ier, &
               limit, lenw, last, iwork, work, za, zb, tx)

   deallocate (iwork, work)
   return
end subroutine dqag_sk3e

subroutine dqag3e(f, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
                  limit, lenw, last, iwork, work, za, zb, tx)
   implicit none
   real*16 a, abserr, b, epsabs, epsrel, f, result, work, tx, za, zb
   integer ier, iwork, key, last, lenw, limit, lvl, l1, l2, l3, neval
   dimension iwork(limit), work(lenw)
   external f
   ier = 6
   neval = 0
   last = 0
   result = 0.0q0
   abserr = 0.0q0
   if (limit .ge. 1 .and. lenw .ge. limit*4) then
      l1 = limit + 1
      l2 = limit + l1
      l3 = limit + l2
      call dqage3e(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
                   ier, work(1), work(l1), work(l2), work(l3), iwork, last, za, zb, tx)
      lvl = 0
   end if
   if (ier .eq. 6) lvl = 1
   if (ier .ne. 0) call xerror("26habnormal return from dqag", 26, ier, lvl)

   return
end subroutine dqag3e

subroutine dqage3e(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                   neval, ier, alist, blist, rlist, elist, iord, last, za, zb, tx)
   implicit none
   real*16 a, abserr, alist, area, area1, area12, area2, a1, a2, b, &
      blist, b1, b2, abs, defabs, defab1, defab2, dmax1, elist, epmach, &
      epsabs, epsrel, errbnd, errmax, error1, error2, erro12, errsum, f, &
      resabs, result, rlist, uflow, za, zb, tx
   integer ier, iord, iroff1, iroff2, k, key, keyf, last, limit, maxerr, neval, &
      nrmax, igt, igk
   dimension alist(limit), blist(limit), elist(limit), iord(limit), &
      rlist(limit)
   external f
   epmach = epsilon(a)
   uflow = tiny(a)
   ier = 0
   neval = 0
   last = 0
   result = 0.0q0
   abserr = 0.0q0
   alist(1) = a
   blist(1) = b
   rlist(1) = 0.0q0
   elist(1) = 0.0q0
   iord(1) = 0
   if (epsabs .le. 0.0q0 .and. &
       epsrel .lt. dmax1(0.5d+02*epmach, 0.5d-28)) ier = 6
   if (ier .ne. 6) then
      keyf = key
      if (key .le. 0) keyf = 1
      if (key .ge. 7) keyf = 6
      neval = 0
      if (keyf .eq. 1) call dqk15_3e(f, a, b, result, abserr, defabs, resabs, za, zb, tx, epsrel)
      if (keyf .eq. 2) call dqk21_3e(f, a, b, result, abserr, defabs, resabs, za, zb, tx, epsrel)
      if (keyf .eq. 3) call dqk31_3e(f, a, b, result, abserr, defabs, resabs, za, zb, tx, epsrel)
      if (keyf .eq. 4) call dqk41_3e(f, a, b, result, abserr, defabs, resabs, za, zb, tx, epsrel)
      if (keyf .eq. 5) call dqk51_3e(f, a, b, result, abserr, defabs, resabs, za, zb, tx, epsrel)
      if (keyf .eq. 6) call dqk61_3e(f, a, b, result, abserr, defabs, resabs, za, zb, tx, epsrel)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      errbnd = dmax1(epsabs, epsrel*abs(result))
      if (abserr .le. 0.5d+02*epmach*defabs .and. abserr .gt. errbnd) ier = 2
      if (limit .eq. 1) ier = 1
      igk = 0
      if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr .ne. resabs) &
          .or. abserr .eq. 0.0q0) igk = 60
      if (igk .ne. 60) then
         errmax = abserr
         maxerr = 1
         area = result
         errsum = abserr
         nrmax = 1
         iroff1 = 0
         iroff2 = 0
         do last = 2, limit
            a1 = alist(maxerr)
            b1 = 0.5q0*(alist(maxerr) + blist(maxerr))
            a2 = b1
            b2 = blist(maxerr)
            if (keyf .eq. 1) call dqk15_3e(f, a1, b1, area1, error1, resabs, defab1, za, zb, tx, epsrel)
            if (keyf .eq. 2) call dqk21_3e(f, a1, b1, area1, error1, resabs, defab1, za, zb, tx, epsrel)
            if (keyf .eq. 3) call dqk31_3e(f, a1, b1, area1, error1, resabs, defab1, za, zb, tx, epsrel)
            if (keyf .eq. 4) call dqk41_3e(f, a1, b1, area1, error1, resabs, defab1, za, zb, tx, epsrel)
            if (keyf .eq. 5) call dqk51_3e(f, a1, b1, area1, error1, resabs, defab1, za, zb, tx, epsrel)
            if (keyf .eq. 6) call dqk61_3e(f, a1, b1, area1, error1, resabs, defab1, za, zb, tx, epsrel)
            if (keyf .eq. 1) call dqk15_3e(f, a2, b2, area2, error2, resabs, defab2, za, zb, tx, epsrel)
            if (keyf .eq. 2) call dqk21_3e(f, a2, b2, area2, error2, resabs, defab2, za, zb, tx, epsrel)
            if (keyf .eq. 3) call dqk31_3e(f, a2, b2, area2, error2, resabs, defab2, za, zb, tx, epsrel)
            if (keyf .eq. 4) call dqk41_3e(f, a2, b2, area2, error2, resabs, defab2, za, zb, tx, epsrel)
            if (keyf .eq. 5) call dqk51_3e(f, a2, b2, area2, error2, resabs, defab2, za, zb, tx, epsrel)
            if (keyf .eq. 6) call dqk61_3e(f, a2, b2, area2, error2, resabs, defab2, za, zb, tx, epsrel)
            neval = neval + 1
            area12 = area1 + area2
            erro12 = error1 + error2
            errsum = errsum + erro12 - errmax
            area = area + area12 - rlist(maxerr)
            if (defab1 .ne. error1 .and. defab2 .ne. error2) then
               if (abs(rlist(maxerr) - area12) .le. 0.1d-04*abs(area12) &
                   .and. erro12 .ge. 0.99q0*errmax) iroff1 = iroff1 + 1
               if (last .gt. 10 .and. erro12 .gt. errmax) iroff2 = iroff2 + 1
            end if
            rlist(maxerr) = area1
            rlist(last) = area2
            errbnd = dmax1(epsabs, epsrel*abs(area))
            if (errsum .gt. errbnd) then
               if (iroff1 .ge. 6 .or. iroff2 .ge. 20) ier = 2
               if (last .eq. limit) ier = 1
               if (dmax1(abs(a1), abs(b2)) .le. (0.1d+01 + 0.1d+03* &
                                                 epmach)*(abs(a2) + 0.1d+04*uflow)) ier = 3
            end if
            igt = 0
            if (error2 .le. error1) then
               alist(last) = a2
               blist(maxerr) = b1
               blist(last) = b2
               elist(maxerr) = error1
               elist(last) = error2
               igt = 20
            end if
            if (igt .ne. 20) then
               alist(maxerr) = a2
               alist(last) = a1
               blist(last) = b1
               rlist(maxerr) = area2
               rlist(last) = area1
               elist(maxerr) = error2
               elist(last) = error1
            end if
            call dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax)
            if (ier .ne. 0 .or. errsum .le. errbnd) exit
         end do
         if (last .eq. limit + 1) last = last - 1
         result = 0.0q0
         do k = 1, last
            result = result + rlist(k)
         end do
         abserr = errsum
      end if
      if (keyf .ne. 1) neval = (10*keyf + 1)*(2*neval + 1)
      if (keyf .eq. 1) neval = 30*neval + 15
   end if
   return
end subroutine dqage3e

subroutine dqk15_3e(f, a, b, result, abserr, resabs, resasc, za, zb, tx, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, za, zb, tx, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(7), fv2(7), wg(4), wgk(8), xgk(8)
   data wg(1)/0.129484966168869693270611432679082q0/
   data wg(2)/0.279705391489276667901467771423780q0/
   data wg(3)/0.381830050505118944950369775488975q0/
   data wg(4)/0.417959183673469387755102040816327q0/

   data xgk(1)/0.991455371120812639206854697526329q0/
   data xgk(2)/0.949107912342758524526189684047851q0/
   data xgk(3)/0.864864423359769072789712788640926q0/
   data xgk(4)/0.741531185599394439863864773280788q0/
   data xgk(5)/0.586087235467691130294144838258730q0/
   data xgk(6)/0.405845151377397166906606412076961q0/
   data xgk(7)/0.207784955007898467600689403773245q0/
   data xgk(8)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.022935322010529224963732008058970q0/
   data wgk(2)/0.063092092629978553290700663189204q0/
   data wgk(3)/0.104790010322250183839876322541518q0/
   data wgk(4)/0.140653259715525918745189590510238q0/
   data wgk(5)/0.169004726639267902826583426598550q0/
   data wgk(6)/0.190350578064785409913256402421014q0/
   data wgk(7)/0.204432940075298892414161999234649q0/
   data wgk(8)/0.209482141084727828012999174891714q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   !fc = f(tx,centr)
   call dqag_sk3f(f, za, zb, eps, fc, ier, 1, centr, tx)
   resg = fc*wg(4)
   resk = fc*wgk(8)
   resabs = abs(resk)
   do j = 1, 3
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 1, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 1, centr + absc, tx)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 4
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 1, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 1, centr + absc, tx)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(8)*abs(fc - reskh)
   do j = 1, 7
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk15_3e

!]]]]]]]]]]]

subroutine dqk21_3e(f, a, b, result, abserr, resabs, resasc, za, zb, tx, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, za, zb, tx, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(10), fv2(10), wg(5), wgk(11), xgk(11)

   data wg(1)/0.066671344308688137593568809893332q0/
   data wg(2)/0.149451349150580593145776339657697q0/
   data wg(3)/0.219086362515982043995534934228163q0/
   data wg(4)/0.269266719309996355091226921569469q0/
   data wg(5)/0.295524224714752870173892994651338q0/

   data xgk(1)/0.995657163025808080735527280689003q0/
   data xgk(2)/0.973906528517171720077964012084452q0/
   data xgk(3)/0.930157491355708226001207180059508q0/
   data xgk(4)/0.865063366688984510732096688423493q0/
   data xgk(5)/0.780817726586416897063717578345042q0/
   data xgk(6)/0.679409568299024406234327365114874q0/
   data xgk(7)/0.562757134668604683339000099272694q0/
   data xgk(8)/0.433395394129247190799265943165784q0/
   data xgk(9)/0.294392862701460198131126603103866q0/
   data xgk(10)/0.148874338981631210884826001129720q0/
   data xgk(11)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.011694638867371874278064396062192q0/
   data wgk(2)/0.032558162307964727478818972459390q0/
   data wgk(3)/0.054755896574351996031381300244580q0/
   data wgk(4)/0.075039674810919952767043140916190q0/
   data wgk(5)/0.093125454583697605535065465083366q0/
   data wgk(6)/0.109387158802297641899210590325805q0/
   data wgk(7)/0.123491976262065851077958109831074q0/
   data wgk(8)/0.134709217311473325928054001771707q0/
   data wgk(9)/0.142775938577060080797094273138717q0/
   data wgk(10)/0.147739104901338491374841515972068q0/
   data wgk(11)/0.149445554002916905664936468389821q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   !fc = f(tx,centr)
   call dqag_sk3f(f, za, zb, eps, fc, ier, 2, centr, tx)

   resk = wgk(11)*fc
   resabs = abs(resk)
   do j = 1, 5
      jtw = 2*j
      absc = hlgth*xgk(jtw)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 2, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 2, centr + absc, tx)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 5
      jtwm1 = 2*j - 1
      absc = hlgth*xgk(jtwm1)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 2, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 2, centr + absc, tx)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(11)*abs(fc - reskh)
   do j = 1, 10
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk21_3e

subroutine dqk31_3e(f, a, b, result, abserr, resabs, resasc, za, zb, tx, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, za, zb, tx, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(15), fv2(15), xgk(16), wgk(16), wg(8)

   data wg(1)/0.030753241996117268354628393577204q0/
   data wg(2)/0.070366047488108124709267416450667q0/
   data wg(3)/0.107159220467171935011869546685869q0/
   data wg(4)/0.139570677926154314447804794511028q0/
   data wg(5)/0.166269205816993933553200860481209q0/
   data wg(6)/0.186161000015562211026800561866423q0/
   data wg(7)/0.198431485327111576456118326443839q0/
   data wg(8)/0.202578241925561272880620199967519q0/

   data xgk(1)/0.998002298693397060285172840152271q0/
   data xgk(2)/0.987992518020485428489565718586613q0/
   data xgk(3)/0.967739075679139134257347978784337q0/
   data xgk(4)/0.937273392400705904307758947710209q0/
   data xgk(5)/0.897264532344081900882509656454496q0/
   data xgk(6)/0.848206583410427216200648320774217q0/
   data xgk(7)/0.790418501442465932967649294817947q0/
   data xgk(8)/0.724417731360170047416186054613938q0/
   data xgk(9)/0.650996741297416970533735895313275q0/
   data xgk(10)/0.570972172608538847537226737253911q0/
   data xgk(11)/0.485081863640239680693655740232351q0/
   data xgk(12)/0.394151347077563369897207370981045q0/
   data xgk(13)/0.299180007153168812166780024266389q0/
   data xgk(14)/0.201194093997434522300628303394596q0/
   data xgk(15)/0.101142066918717499027074231447392q0/
   data xgk(16)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.005377479872923348987792051430128q0/
   data wgk(2)/0.015007947329316122538374763075807q0/
   data wgk(3)/0.025460847326715320186874001019653q0/
   data wgk(4)/0.035346360791375846222037948478360q0/
   data wgk(5)/0.044589751324764876608227299373280q0/
   data wgk(6)/0.053481524690928087265343147239430q0/
   data wgk(7)/0.062009567800670640285139230960803q0/
   data wgk(8)/0.069854121318728258709520077099147q0/
   data wgk(9)/0.076849680757720378894432777482659q0/
   data wgk(10)/0.083080502823133021038289247286104q0/
   data wgk(11)/0.088564443056211770647275443693774q0/
   data wgk(12)/0.093126598170825321225486872747346q0/
   data wgk(13)/0.096642726983623678505179907627589q0/
   data wgk(14)/0.099173598721791959332393173484603q0/
   data wgk(15)/0.100769845523875595044946662617570q0/
   data wgk(16)/0.101330007014791549017374792767493q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)
   !fc = f(tx,centr)
   call dqag_sk3f(f, za, zb, eps, fc, ier, 3, centr, tx)

   resg = wg(8)*fc
   resk = wgk(16)*fc
   resabs = abs(resk)
   do j = 1, 7
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 3, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 3, centr + absc, tx)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 8
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 3, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 3, centr + absc, tx)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(16)*abs(fc - reskh)
   do j = 1, 15
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk31_3e

subroutine dqk41_3e(f, a, b, result, abserr, resabs, resasc, za, zb, tx, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, za, zb, tx, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(20), fv2(20), xgk(21), wgk(21), wg(10)

   data wg(1)/0.017614007139152118311861962351853q0/
   data wg(2)/0.040601429800386941331039952274932q0/
   data wg(3)/0.062672048334109063569506535187042q0/
   data wg(4)/0.083276741576704748724758143222046q0/
   data wg(5)/0.101930119817240435036750135480350q0/
   data wg(6)/0.118194531961518417312377377711382q0/
   data wg(7)/0.131688638449176626898494499748163q0/
   data wg(8)/0.142096109318382051329298325067165q0/
   data wg(9)/0.149172986472603746787828737001969q0/
   data wg(10)/0.152753387130725850698084331955098q0/
   !
   data xgk(1)/0.998859031588277663838315576545863q0/
   data xgk(2)/0.993128599185094924786122388471320q0/
   data xgk(3)/0.981507877450250259193342994720217q0/
   data xgk(4)/0.963971927277913791267666131197277q0/
   data xgk(5)/0.940822633831754753519982722212443q0/
   data xgk(6)/0.912234428251325905867752441203298q0/
   data xgk(7)/0.878276811252281976077442995113078q0/
   data xgk(8)/0.839116971822218823394529061701521q0/
   data xgk(9)/0.795041428837551198350638833272788q0/
   data xgk(10)/0.746331906460150792614305070355642q0/
   data xgk(11)/0.693237656334751384805490711845932q0/
   data xgk(12)/0.636053680726515025452836696226286q0/
   data xgk(13)/0.575140446819710315342946036586425q0/
   data xgk(14)/0.510867001950827098004364050955251q0/
   data xgk(15)/0.443593175238725103199992213492640q0/
   data xgk(16)/0.373706088715419560672548177024927q0/
   data xgk(17)/0.301627868114913004320555356858592q0/
   data xgk(18)/0.227785851141645078080496195368575q0/
   data xgk(19)/0.152605465240922675505220241022678q0/
   data xgk(20)/0.076526521133497333754640409398838q0/
   data xgk(21)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.003073583718520531501218293246031q0/
   data wgk(2)/0.008600269855642942198661787950102q0/
   data wgk(3)/0.014626169256971252983787960308868q0/
   data wgk(4)/0.020388373461266523598010231432755q0/
   data wgk(5)/0.025882133604951158834505067096153q0/
   data wgk(6)/0.031287306777032798958543119323801q0/
   data wgk(7)/0.036600169758200798030557240707211q0/
   data wgk(8)/0.041668873327973686263788305936895q0/
   data wgk(9)/0.046434821867497674720231880926108q0/
   data wgk(10)/0.050944573923728691932707670050345q0/
   data wgk(11)/0.055195105348285994744832372419777q0/
   data wgk(12)/0.059111400880639572374967220648594q0/
   data wgk(13)/0.062653237554781168025870122174255q0/
   data wgk(14)/0.065834597133618422111563556969398q0/
   data wgk(15)/0.068648672928521619345623411885368q0/
   data wgk(16)/0.071054423553444068305790361723210q0/
   data wgk(17)/0.073030690332786667495189417658913q0/
   data wgk(18)/0.074582875400499188986581418362488q0/
   data wgk(19)/0.075704497684556674659542775376617q0/
   data wgk(20)/0.076377867672080736705502835038061q0/
   data wgk(21)/0.076600711917999656445049901530102q0/

   epmach = epsilon(a)
   uflow = tiny(a)

   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   !fc = f(tx,centr)
   call dqag_sk3f(f, za, zb, eps, fc, ier, 4, centr, tx)

   resk = wgk(21)*fc
   resabs = abs(resk)
   do j = 1, 10
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 4, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 4, centr + absc, tx)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 10
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 4, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 4, centr + absc, tx)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(21)*abs(fc - reskh)
   do j = 1, 20
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk41_3e

subroutine dqk51_3e(f, a, b, result, abserr, resabs, resasc, za, zb, tx, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, za, zb, tx, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(25), fv2(25), xgk(26), wgk(26), wg(13)

   data wg(1)/0.011393798501026287947902964113235q0/
   data wg(2)/0.026354986615032137261901815295299q0/
   data wg(3)/0.040939156701306312655623487711646q0/
   data wg(4)/0.054904695975835191925936891540473q0/
   data wg(5)/0.068038333812356917207187185656708q0/
   data wg(6)/0.080140700335001018013234959669111q0/
   data wg(7)/0.091028261982963649811497220702892q0/
   data wg(8)/0.100535949067050644202206890392686q0/
   data wg(9)/0.108519624474263653116093957050117q0/
   data wg(10)/0.114858259145711648339325545869556q0/
   data wg(11)/0.119455763535784772228178126512901q0/
   data wg(12)/0.122242442990310041688959518945852q0/
   data wg(13)/0.123176053726715451203902873079050q0/
   !
   data xgk(1)/0.999262104992609834193457486540341q0/
   data xgk(2)/0.995556969790498097908784946893902q0/
   data xgk(3)/0.988035794534077247637331014577406q0/
   data xgk(4)/0.976663921459517511498315386479594q0/
   data xgk(5)/0.961614986425842512418130033660167q0/
   data xgk(6)/0.942974571228974339414011169658471q0/
   data xgk(7)/0.920747115281701561746346084546331q0/
   data xgk(8)/0.894991997878275368851042006782805q0/
   data xgk(9)/0.865847065293275595448996969588340q0/
   data xgk(10)/0.833442628760834001421021108693570q0/
   data xgk(11)/0.797873797998500059410410904994307q0/
   data xgk(12)/0.759259263037357630577282865204361q0/
   data xgk(13)/0.717766406813084388186654079773298q0/
   data xgk(14)/0.673566368473468364485120633247622q0/
   data xgk(15)/0.626810099010317412788122681624518q0/
   data xgk(16)/0.577662930241222967723689841612654q0/
   data xgk(17)/0.526325284334719182599623778158010q0/
   data xgk(18)/0.473002731445714960522182115009192q0/
   data xgk(19)/0.417885382193037748851814394594572q0/
   data xgk(20)/0.361172305809387837735821730127641q0/
   data xgk(21)/0.303089538931107830167478909980339q0/
   data xgk(22)/0.243866883720988432045190362797452q0/
   data xgk(23)/0.183718939421048892015969888759528q0/
   data xgk(24)/0.122864692610710396387359818808037q0/
   data xgk(25)/0.061544483005685078886546392366797q0/
   data xgk(26)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.001987383892330315926507851882843q0/
   data wgk(2)/0.005561932135356713758040236901066q0/
   data wgk(3)/0.009473973386174151607207710523655q0/
   data wgk(4)/0.013236229195571674813656405846976q0/
   data wgk(5)/0.016847817709128298231516667536336q0/
   data wgk(6)/0.020435371145882835456568292235939q0/
   data wgk(7)/0.024009945606953216220092489164881q0/
   data wgk(8)/0.027475317587851737802948455517811q0/
   data wgk(9)/0.030792300167387488891109020215229q0/
   data wgk(10)/0.034002130274329337836748795229551q0/
   data wgk(11)/0.037116271483415543560330625367620q0/
   data wgk(12)/0.040083825504032382074839284467076q0/
   data wgk(13)/0.042872845020170049476895792439495q0/
   data wgk(14)/0.045502913049921788909870584752660q0/
   data wgk(15)/0.047982537138836713906392255756915q0/
   data wgk(16)/0.050277679080715671963325259433440q0/
   data wgk(17)/0.052362885806407475864366712137873q0/
   data wgk(18)/0.054251129888545490144543370459876q0/
   data wgk(19)/0.055950811220412317308240686382747q0/
   data wgk(20)/0.057437116361567832853582693939506q0/
   data wgk(21)/0.058689680022394207961974175856788q0/
   data wgk(22)/0.059720340324174059979099291932562q0/
   data wgk(23)/0.060539455376045862945360267517565q0/
   data wgk(24)/0.061128509717053048305859030416293q0/
   data wgk(25)/0.061471189871425316661544131965264q0/
   !       note: wgk (26) was calculated from the values of wgk(1..25)
   data wgk(26)/0.061580818067832935078759824240066q0/

   epmach = epsilon(a)
   uflow = tiny(a)

   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)
   !fc = f(tx,centr)
   call dqag_sk3f(f, za, zb, eps, fc, ier, 5, centr, tx)

   resg = wg(13)*fc
   resk = wgk(26)*fc
   resabs = abs(resk)
   do j = 1, 12
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 5, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 5, centr + absc, tx)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 13
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 5, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 5, centr + absc, tx)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(26)*abs(fc - reskh)
   do j = 1, 25
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk51_3e

subroutine dqk61_3e(f, a, b, result, abserr, resabs, resasc, za, zb, tx, eps)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, za, zb, tx, eps
   integer j, jtw, jtwm1, ier
   external f
   dimension fv1(30), fv2(30), xgk(31), wgk(31), wg(15)

   data wg(1)/0.007968192496166605615465883474674q0/
   data wg(2)/0.018466468311090959142302131912047q0/
   data wg(3)/0.028784707883323369349719179611292q0/
   data wg(4)/0.038799192569627049596801936446348q0/
   data wg(5)/0.048402672830594052902938140422808q0/
   data wg(6)/0.057493156217619066481721689402056q0/
   data wg(7)/0.065974229882180495128128515115962q0/
   data wg(8)/0.073755974737705206268243850022191q0/
   data wg(9)/0.080755895229420215354694938460530q0/
   data wg(10)/0.086899787201082979802387530715126q0/
   data wg(11)/0.092122522237786128717632707087619q0/
   data wg(12)/0.096368737174644259639468626351810q0/
   data wg(13)/0.099593420586795267062780282103569q0/
   data wg(14)/0.101762389748405504596428952168554q0/
   data wg(15)/0.102852652893558840341285636705415q0/
   !
   data xgk(1)/0.999484410050490637571325895705811q0/
   data xgk(2)/0.996893484074649540271630050918695q0/
   data xgk(3)/0.991630996870404594858628366109486q0/
   data xgk(4)/0.983668123279747209970032581605663q0/
   data xgk(5)/0.973116322501126268374693868423707q0/
   data xgk(6)/0.960021864968307512216871025581798q0/
   data xgk(7)/0.944374444748559979415831324037439q0/
   data xgk(8)/0.926200047429274325879324277080474q0/
   data xgk(9)/0.905573307699907798546522558925958q0/
   data xgk(10)/0.882560535792052681543116462530226q0/
   data xgk(11)/0.857205233546061098958658510658944q0/
   data xgk(12)/0.829565762382768397442898119732502q0/
   data xgk(13)/0.799727835821839083013668942322683q0/
   data xgk(14)/0.767777432104826194917977340974503q0/
   data xgk(15)/0.733790062453226804726171131369528q0/
   data xgk(16)/0.697850494793315796932292388026640q0/
   data xgk(17)/0.660061064126626961370053668149271q0/
   data xgk(18)/0.620526182989242861140477556431189q0/
   data xgk(19)/0.579345235826361691756024932172540q0/
   data xgk(20)/0.536624148142019899264169793311073q0/
   data xgk(21)/0.492480467861778574993693061207709q0/
   data xgk(22)/0.447033769538089176780609900322854q0/
   data xgk(23)/0.400401254830394392535476211542661q0/
   data xgk(24)/0.352704725530878113471037207089374q0/
   data xgk(25)/0.304073202273625077372677107199257q0/
   data xgk(26)/0.254636926167889846439805129817805q0/
   data xgk(27)/0.204525116682309891438957671002025q0/
   data xgk(28)/0.153869913608583546963794672743256q0/
   data xgk(29)/0.102806937966737030147096751318001q0/
   data xgk(30)/0.051471842555317695833025213166723q0/
   data xgk(31)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.001389013698677007624551591226760q0/
   data wgk(2)/0.003890461127099884051267201844516q0/
   data wgk(3)/0.006630703915931292173319826369750q0/
   data wgk(4)/0.009273279659517763428441146892024q0/
   data wgk(5)/0.011823015253496341742232898853251q0/
   data wgk(6)/0.014369729507045804812451432443580q0/
   data wgk(7)/0.016920889189053272627572289420322q0/
   data wgk(8)/0.019414141193942381173408951050128q0/
   data wgk(9)/0.021828035821609192297167485738339q0/
   data wgk(10)/0.024191162078080601365686370725232q0/
   data wgk(11)/0.026509954882333101610601709335075q0/
   data wgk(12)/0.028754048765041292843978785354334q0/
   data wgk(13)/0.030907257562387762472884252943092q0/
   data wgk(14)/0.032981447057483726031814191016854q0/
   data wgk(15)/0.034979338028060024137499670731468q0/
   data wgk(16)/0.036882364651821229223911065617136q0/
   data wgk(17)/0.038678945624727592950348651532281q0/
   data wgk(18)/0.040374538951535959111995279752468q0/
   data wgk(19)/0.041969810215164246147147541285970q0/
   data wgk(20)/0.043452539701356069316831728117073q0/
   data wgk(21)/0.044814800133162663192355551616723q0/
   data wgk(22)/0.046059238271006988116271735559374q0/
   data wgk(23)/0.047185546569299153945261478181099q0/
   data wgk(24)/0.048185861757087129140779492298305q0/
   data wgk(25)/0.049055434555029778887528165367238q0/
   data wgk(26)/0.049795683427074206357811569379942q0/
   data wgk(27)/0.050405921402782346840893085653585q0/
   data wgk(28)/0.050881795898749606492297473049805q0/
   data wgk(29)/0.051221547849258772170656282604944q0/
   data wgk(30)/0.051426128537459025933862879215781q0/
   data wgk(31)/0.051494729429451567558340433647099q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   !
   centr = 0.5q0*(b + a)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   !fc = f(tx,centr)
   call dqag_sk3f(f, za, zb, eps, fc, ier, 6, centr, tx)
   resk = wgk(31)*fc
   resabs = abs(resk)
   do j = 1, 15
      jtw = j*2
      absc = hlgth*xgk(jtw)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 6, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 6, centr + absc, tx)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 15
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      !fval1 = f(tx,centr-absc)
      !fval2 = f(tx,centr+absc)
      call dqag_sk3f(f, za, zb, eps, fval1, ier, 6, centr - absc, tx)
      call dqag_sk3f(f, za, zb, eps, fval2, ier, 6, centr + absc, tx)

      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(31)*abs(fc - reskh)
   do j = 1, 30
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk61_3e

!]]]]]]]]]]

subroutine dqag_sk3f(f, a, b, eps, s, ier, key, ty, tx)
   implicit none
   integer, intent(in)::key
   real*16, intent(in)::a, b, eps, tx, ty
   real*16, intent(out)::s
   integer, intent(out)::ier
   real*16::f
   ! ier = 0 : success, converged.
   ! ier > 0 : fail, not converged.

   integer::neval, limit, lenw, last
   integer, allocatable::iwork(:)
   real*16, allocatable::work(:)
   real*16::epsabs, abserr
   external::f

   ier = 1
   epsabs = -1q0
   !limit=200
   limit = 1000
   lenw = limit*4
   allocate (iwork(1:limit), work(1:lenw))

   s = 0q0
   call dqag3f(f, a, b, epsabs, eps, key, s, abserr, neval, ier, &
               limit, lenw, last, iwork, work, ty, tx)

   deallocate (iwork, work)
   return
end subroutine dqag_sk3f

subroutine dqag3f(f, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
                  limit, lenw, last, iwork, work, ty, tx)
   implicit none
   real*16 a, abserr, b, epsabs, epsrel, f, result, work, tx, ty
   integer ier, iwork, key, last, lenw, limit, lvl, l1, l2, l3, neval
   dimension iwork(limit), work(lenw)
   external f
   ier = 6
   neval = 0
   last = 0
   result = 0.0q0
   abserr = 0.0q0
   if (limit .ge. 1 .and. lenw .ge. limit*4) then
      l1 = limit + 1
      l2 = limit + l1
      l3 = limit + l2
      call dqage3f(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
                   ier, work(1), work(l1), work(l2), work(l3), iwork, last, ty, tx)
      lvl = 0
   end if
   if (ier .eq. 6) lvl = 1
   if (ier .ne. 0) call xerror("26habnormal return from dqag", 26, ier, lvl)

   return
end subroutine dqag3f

subroutine dqage3f(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                   neval, ier, alist, blist, rlist, elist, iord, last, ty, tx)
   implicit none
   real*16 a, abserr, alist, area, area1, area12, area2, a1, a2, b, &
      blist, b1, b2, abs, defabs, defab1, defab2, dmax1, elist, epmach, &
      epsabs, epsrel, errbnd, errmax, error1, error2, erro12, errsum, f, &
      resabs, result, rlist, uflow, ty, tx
   integer ier, iord, iroff1, iroff2, k, key, keyf, last, limit, maxerr, neval, &
      nrmax, igt, igk
   dimension alist(limit), blist(limit), elist(limit), iord(limit), &
      rlist(limit)
   external f
   epmach = epsilon(a)
   uflow = tiny(a)
   ier = 0
   neval = 0
   last = 0
   result = 0.0q0
   abserr = 0.0q0
   alist(1) = a
   blist(1) = b
   rlist(1) = 0.0q0
   elist(1) = 0.0q0
   iord(1) = 0
   if (epsabs .le. 0.0q0 .and. &
       epsrel .lt. dmax1(0.5d+02*epmach, 0.5d-28)) ier = 6
   if (ier .ne. 6) then
      keyf = key
      if (key .le. 0) keyf = 1
      if (key .ge. 7) keyf = 6
      neval = 0
      if (keyf .eq. 1) call dqk15_3f(f, a, b, result, abserr, defabs, resabs, ty, tx)
      if (keyf .eq. 2) call dqk21_3f(f, a, b, result, abserr, defabs, resabs, ty, tx)
      if (keyf .eq. 3) call dqk31_3f(f, a, b, result, abserr, defabs, resabs, ty, tx)
      if (keyf .eq. 4) call dqk41_3f(f, a, b, result, abserr, defabs, resabs, ty, tx)
      if (keyf .eq. 5) call dqk51_3f(f, a, b, result, abserr, defabs, resabs, ty, tx)
      if (keyf .eq. 6) call dqk61_3f(f, a, b, result, abserr, defabs, resabs, ty, tx)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      errbnd = dmax1(epsabs, epsrel*abs(result))
      if (abserr .le. 0.5d+02*epmach*defabs .and. abserr .gt. errbnd) ier = 2
      if (limit .eq. 1) ier = 1
      igk = 0
      if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr .ne. resabs) &
          .or. abserr .eq. 0.0q0) igk = 60
      if (igk .ne. 60) then
         errmax = abserr
         maxerr = 1
         area = result
         errsum = abserr
         nrmax = 1
         iroff1 = 0
         iroff2 = 0
         do last = 2, limit
            a1 = alist(maxerr)
            b1 = 0.5q0*(alist(maxerr) + blist(maxerr))
            a2 = b1
            b2 = blist(maxerr)
            if (keyf .eq. 1) call dqk15_3f(f, a1, b1, area1, error1, resabs, defab1, ty, tx)
            if (keyf .eq. 2) call dqk21_3f(f, a1, b1, area1, error1, resabs, defab1, ty, tx)
            if (keyf .eq. 3) call dqk31_3f(f, a1, b1, area1, error1, resabs, defab1, ty, tx)
            if (keyf .eq. 4) call dqk41_3f(f, a1, b1, area1, error1, resabs, defab1, ty, tx)
            if (keyf .eq. 5) call dqk51_3f(f, a1, b1, area1, error1, resabs, defab1, ty, tx)
            if (keyf .eq. 6) call dqk61_3f(f, a1, b1, area1, error1, resabs, defab1, ty, tx)
            if (keyf .eq. 1) call dqk15_3f(f, a2, b2, area2, error2, resabs, defab2, ty, tx)
            if (keyf .eq. 2) call dqk21_3f(f, a2, b2, area2, error2, resabs, defab2, ty, tx)
            if (keyf .eq. 3) call dqk31_3f(f, a2, b2, area2, error2, resabs, defab2, ty, tx)
            if (keyf .eq. 4) call dqk41_3f(f, a2, b2, area2, error2, resabs, defab2, ty, tx)
            if (keyf .eq. 5) call dqk51_3f(f, a2, b2, area2, error2, resabs, defab2, ty, tx)
            if (keyf .eq. 6) call dqk61_3f(f, a2, b2, area2, error2, resabs, defab2, ty, tx)
            neval = neval + 1
            area12 = area1 + area2
            erro12 = error1 + error2
            errsum = errsum + erro12 - errmax
            area = area + area12 - rlist(maxerr)
            if (defab1 .ne. error1 .and. defab2 .ne. error2) then
               if (abs(rlist(maxerr) - area12) .le. 0.1d-04*abs(area12) &
                   .and. erro12 .ge. 0.99q0*errmax) iroff1 = iroff1 + 1
               if (last .gt. 10 .and. erro12 .gt. errmax) iroff2 = iroff2 + 1
            end if
            rlist(maxerr) = area1
            rlist(last) = area2
            errbnd = dmax1(epsabs, epsrel*abs(area))
            if (errsum .gt. errbnd) then
               if (iroff1 .ge. 6 .or. iroff2 .ge. 20) ier = 2
               if (last .eq. limit) ier = 1
               if (dmax1(abs(a1), abs(b2)) .le. (0.1d+01 + 0.1d+03* &
                                                 epmach)*(abs(a2) + 0.1d+04*uflow)) ier = 3
            end if
            igt = 0
            if (error2 .le. error1) then
               alist(last) = a2
               blist(maxerr) = b1
               blist(last) = b2
               elist(maxerr) = error1
               elist(last) = error2
               igt = 20
            end if
            if (igt .ne. 20) then
               alist(maxerr) = a2
               alist(last) = a1
               blist(last) = b1
               rlist(maxerr) = area2
               rlist(last) = area1
               elist(maxerr) = error2
               elist(last) = error1
            end if
            call dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax)
            if (ier .ne. 0 .or. errsum .le. errbnd) exit
         end do
         if (last .eq. limit + 1) last = last - 1
         result = 0.0q0
         do k = 1, last
            result = result + rlist(k)
         end do
         abserr = errsum
      end if
      if (keyf .ne. 1) neval = (10*keyf + 1)*(2*neval + 1)
      if (keyf .eq. 1) neval = 30*neval + 15
   end if
   return
end subroutine dqage3f

subroutine dqk15_3f(f, a, b, result, abserr, resabs, resasc, ty, tx)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, ty, tx
   integer j, jtw, jtwm1
   external f
   dimension fv1(7), fv2(7), wg(4), wgk(8), xgk(8)
   data wg(1)/0.129484966168869693270611432679082q0/
   data wg(2)/0.279705391489276667901467771423780q0/
   data wg(3)/0.381830050505118944950369775488975q0/
   data wg(4)/0.417959183673469387755102040816327q0/

   data xgk(1)/0.991455371120812639206854697526329q0/
   data xgk(2)/0.949107912342758524526189684047851q0/
   data xgk(3)/0.864864423359769072789712788640926q0/
   data xgk(4)/0.741531185599394439863864773280788q0/
   data xgk(5)/0.586087235467691130294144838258730q0/
   data xgk(6)/0.405845151377397166906606412076961q0/
   data xgk(7)/0.207784955007898467600689403773245q0/
   data xgk(8)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.022935322010529224963732008058970q0/
   data wgk(2)/0.063092092629978553290700663189204q0/
   data wgk(3)/0.104790010322250183839876322541518q0/
   data wgk(4)/0.140653259715525918745189590510238q0/
   data wgk(5)/0.169004726639267902826583426598550q0/
   data wgk(6)/0.190350578064785409913256402421014q0/
   data wgk(7)/0.204432940075298892414161999234649q0/
   data wgk(8)/0.209482141084727828012999174891714q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   fc = f(tx, ty, centr)
   resg = fc*wg(4)
   resk = fc*wgk(8)
   resabs = abs(resk)
   do j = 1, 3
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)

      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 4
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(8)*abs(fc - reskh)
   do j = 1, 7
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk15_3f

!######

subroutine dqk21_3f(f, a, b, result, abserr, resabs, resasc, ty, tx)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, tx, ty
   integer j, jtw, jtwm1
   external f
   dimension fv1(10), fv2(10), wg(5), wgk(11), xgk(11)

   data wg(1)/0.066671344308688137593568809893332q0/
   data wg(2)/0.149451349150580593145776339657697q0/
   data wg(3)/0.219086362515982043995534934228163q0/
   data wg(4)/0.269266719309996355091226921569469q0/
   data wg(5)/0.295524224714752870173892994651338q0/

   data xgk(1)/0.995657163025808080735527280689003q0/
   data xgk(2)/0.973906528517171720077964012084452q0/
   data xgk(3)/0.930157491355708226001207180059508q0/
   data xgk(4)/0.865063366688984510732096688423493q0/
   data xgk(5)/0.780817726586416897063717578345042q0/
   data xgk(6)/0.679409568299024406234327365114874q0/
   data xgk(7)/0.562757134668604683339000099272694q0/
   data xgk(8)/0.433395394129247190799265943165784q0/
   data xgk(9)/0.294392862701460198131126603103866q0/
   data xgk(10)/0.148874338981631210884826001129720q0/
   data xgk(11)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.011694638867371874278064396062192q0/
   data wgk(2)/0.032558162307964727478818972459390q0/
   data wgk(3)/0.054755896574351996031381300244580q0/
   data wgk(4)/0.075039674810919952767043140916190q0/
   data wgk(5)/0.093125454583697605535065465083366q0/
   data wgk(6)/0.109387158802297641899210590325805q0/
   data wgk(7)/0.123491976262065851077958109831074q0/
   data wgk(8)/0.134709217311473325928054001771707q0/
   data wgk(9)/0.142775938577060080797094273138717q0/
   data wgk(10)/0.147739104901338491374841515972068q0/
   data wgk(11)/0.149445554002916905664936468389821q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   fc = f(tx, ty, centr)
   resk = wgk(11)*fc
   resabs = abs(resk)
   do j = 1, 5
      jtw = 2*j
      absc = hlgth*xgk(jtw)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 5
      jtwm1 = 2*j - 1
      absc = hlgth*xgk(jtwm1)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(11)*abs(fc - reskh)
   do j = 1, 10
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk21_3f

subroutine dqk31_3f(f, a, b, result, abserr, resabs, resasc, ty, tx)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, tx, ty
   integer j, jtw, jtwm1
   external f
   dimension fv1(15), fv2(15), xgk(16), wgk(16), wg(8)

   data wg(1)/0.030753241996117268354628393577204q0/
   data wg(2)/0.070366047488108124709267416450667q0/
   data wg(3)/0.107159220467171935011869546685869q0/
   data wg(4)/0.139570677926154314447804794511028q0/
   data wg(5)/0.166269205816993933553200860481209q0/
   data wg(6)/0.186161000015562211026800561866423q0/
   data wg(7)/0.198431485327111576456118326443839q0/
   data wg(8)/0.202578241925561272880620199967519q0/

   data xgk(1)/0.998002298693397060285172840152271q0/
   data xgk(2)/0.987992518020485428489565718586613q0/
   data xgk(3)/0.967739075679139134257347978784337q0/
   data xgk(4)/0.937273392400705904307758947710209q0/
   data xgk(5)/0.897264532344081900882509656454496q0/
   data xgk(6)/0.848206583410427216200648320774217q0/
   data xgk(7)/0.790418501442465932967649294817947q0/
   data xgk(8)/0.724417731360170047416186054613938q0/
   data xgk(9)/0.650996741297416970533735895313275q0/
   data xgk(10)/0.570972172608538847537226737253911q0/
   data xgk(11)/0.485081863640239680693655740232351q0/
   data xgk(12)/0.394151347077563369897207370981045q0/
   data xgk(13)/0.299180007153168812166780024266389q0/
   data xgk(14)/0.201194093997434522300628303394596q0/
   data xgk(15)/0.101142066918717499027074231447392q0/
   data xgk(16)/0.000000000000000000000000000000000q0/

   data wgk(1)/0.005377479872923348987792051430128q0/
   data wgk(2)/0.015007947329316122538374763075807q0/
   data wgk(3)/0.025460847326715320186874001019653q0/
   data wgk(4)/0.035346360791375846222037948478360q0/
   data wgk(5)/0.044589751324764876608227299373280q0/
   data wgk(6)/0.053481524690928087265343147239430q0/
   data wgk(7)/0.062009567800670640285139230960803q0/
   data wgk(8)/0.069854121318728258709520077099147q0/
   data wgk(9)/0.076849680757720378894432777482659q0/
   data wgk(10)/0.083080502823133021038289247286104q0/
   data wgk(11)/0.088564443056211770647275443693774q0/
   data wgk(12)/0.093126598170825321225486872747346q0/
   data wgk(13)/0.096642726983623678505179907627589q0/
   data wgk(14)/0.099173598721791959332393173484603q0/
   data wgk(15)/0.100769845523875595044946662617570q0/
   data wgk(16)/0.101330007014791549017374792767493q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)
   fc = f(tx, ty, centr)
   resg = wg(8)*fc
   resk = wgk(16)*fc
   resabs = abs(resk)
   do j = 1, 7
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 8
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(16)*abs(fc - reskh)
   do j = 1, 15
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk31_3f

subroutine dqk41_3f(f, a, b, result, abserr, resabs, resasc, ty, tx)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, tx, ty
   integer j, jtw, jtwm1
   external f
   dimension fv1(20), fv2(20), xgk(21), wgk(21), wg(10)

   data wg(1)/0.017614007139152118311861962351853q0/
   data wg(2)/0.040601429800386941331039952274932q0/
   data wg(3)/0.062672048334109063569506535187042q0/
   data wg(4)/0.083276741576704748724758143222046q0/
   data wg(5)/0.101930119817240435036750135480350q0/
   data wg(6)/0.118194531961518417312377377711382q0/
   data wg(7)/0.131688638449176626898494499748163q0/
   data wg(8)/0.142096109318382051329298325067165q0/
   data wg(9)/0.149172986472603746787828737001969q0/
   data wg(10)/0.152753387130725850698084331955098q0/
   !
   data xgk(1)/0.998859031588277663838315576545863q0/
   data xgk(2)/0.993128599185094924786122388471320q0/
   data xgk(3)/0.981507877450250259193342994720217q0/
   data xgk(4)/0.963971927277913791267666131197277q0/
   data xgk(5)/0.940822633831754753519982722212443q0/
   data xgk(6)/0.912234428251325905867752441203298q0/
   data xgk(7)/0.878276811252281976077442995113078q0/
   data xgk(8)/0.839116971822218823394529061701521q0/
   data xgk(9)/0.795041428837551198350638833272788q0/
   data xgk(10)/0.746331906460150792614305070355642q0/
   data xgk(11)/0.693237656334751384805490711845932q0/
   data xgk(12)/0.636053680726515025452836696226286q0/
   data xgk(13)/0.575140446819710315342946036586425q0/
   data xgk(14)/0.510867001950827098004364050955251q0/
   data xgk(15)/0.443593175238725103199992213492640q0/
   data xgk(16)/0.373706088715419560672548177024927q0/
   data xgk(17)/0.301627868114913004320555356858592q0/
   data xgk(18)/0.227785851141645078080496195368575q0/
   data xgk(19)/0.152605465240922675505220241022678q0/
   data xgk(20)/0.076526521133497333754640409398838q0/
   data xgk(21)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.003073583718520531501218293246031q0/
   data wgk(2)/0.008600269855642942198661787950102q0/
   data wgk(3)/0.014626169256971252983787960308868q0/
   data wgk(4)/0.020388373461266523598010231432755q0/
   data wgk(5)/0.025882133604951158834505067096153q0/
   data wgk(6)/0.031287306777032798958543119323801q0/
   data wgk(7)/0.036600169758200798030557240707211q0/
   data wgk(8)/0.041668873327973686263788305936895q0/
   data wgk(9)/0.046434821867497674720231880926108q0/
   data wgk(10)/0.050944573923728691932707670050345q0/
   data wgk(11)/0.055195105348285994744832372419777q0/
   data wgk(12)/0.059111400880639572374967220648594q0/
   data wgk(13)/0.062653237554781168025870122174255q0/
   data wgk(14)/0.065834597133618422111563556969398q0/
   data wgk(15)/0.068648672928521619345623411885368q0/
   data wgk(16)/0.071054423553444068305790361723210q0/
   data wgk(17)/0.073030690332786667495189417658913q0/
   data wgk(18)/0.074582875400499188986581418362488q0/
   data wgk(19)/0.075704497684556674659542775376617q0/
   data wgk(20)/0.076377867672080736705502835038061q0/
   data wgk(21)/0.076600711917999656445049901530102q0/

   epmach = epsilon(a)
   uflow = tiny(a)

   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   fc = f(tx, ty, centr)
   resk = wgk(21)*fc
   resabs = abs(resk)
   do j = 1, 10
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 10
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(21)*abs(fc - reskh)
   do j = 1, 20
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk41_3f

subroutine dqk51_3f(f, a, b, result, abserr, resabs, resasc, ty, tx)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, tx, ty
   integer j, jtw, jtwm1
   external f
   dimension fv1(25), fv2(25), xgk(26), wgk(26), wg(13)

   data wg(1)/0.011393798501026287947902964113235q0/
   data wg(2)/0.026354986615032137261901815295299q0/
   data wg(3)/0.040939156701306312655623487711646q0/
   data wg(4)/0.054904695975835191925936891540473q0/
   data wg(5)/0.068038333812356917207187185656708q0/
   data wg(6)/0.080140700335001018013234959669111q0/
   data wg(7)/0.091028261982963649811497220702892q0/
   data wg(8)/0.100535949067050644202206890392686q0/
   data wg(9)/0.108519624474263653116093957050117q0/
   data wg(10)/0.114858259145711648339325545869556q0/
   data wg(11)/0.119455763535784772228178126512901q0/
   data wg(12)/0.122242442990310041688959518945852q0/
   data wg(13)/0.123176053726715451203902873079050q0/
   !
   data xgk(1)/0.999262104992609834193457486540341q0/
   data xgk(2)/0.995556969790498097908784946893902q0/
   data xgk(3)/0.988035794534077247637331014577406q0/
   data xgk(4)/0.976663921459517511498315386479594q0/
   data xgk(5)/0.961614986425842512418130033660167q0/
   data xgk(6)/0.942974571228974339414011169658471q0/
   data xgk(7)/0.920747115281701561746346084546331q0/
   data xgk(8)/0.894991997878275368851042006782805q0/
   data xgk(9)/0.865847065293275595448996969588340q0/
   data xgk(10)/0.833442628760834001421021108693570q0/
   data xgk(11)/0.797873797998500059410410904994307q0/
   data xgk(12)/0.759259263037357630577282865204361q0/
   data xgk(13)/0.717766406813084388186654079773298q0/
   data xgk(14)/0.673566368473468364485120633247622q0/
   data xgk(15)/0.626810099010317412788122681624518q0/
   data xgk(16)/0.577662930241222967723689841612654q0/
   data xgk(17)/0.526325284334719182599623778158010q0/
   data xgk(18)/0.473002731445714960522182115009192q0/
   data xgk(19)/0.417885382193037748851814394594572q0/
   data xgk(20)/0.361172305809387837735821730127641q0/
   data xgk(21)/0.303089538931107830167478909980339q0/
   data xgk(22)/0.243866883720988432045190362797452q0/
   data xgk(23)/0.183718939421048892015969888759528q0/
   data xgk(24)/0.122864692610710396387359818808037q0/
   data xgk(25)/0.061544483005685078886546392366797q0/
   data xgk(26)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.001987383892330315926507851882843q0/
   data wgk(2)/0.005561932135356713758040236901066q0/
   data wgk(3)/0.009473973386174151607207710523655q0/
   data wgk(4)/0.013236229195571674813656405846976q0/
   data wgk(5)/0.016847817709128298231516667536336q0/
   data wgk(6)/0.020435371145882835456568292235939q0/
   data wgk(7)/0.024009945606953216220092489164881q0/
   data wgk(8)/0.027475317587851737802948455517811q0/
   data wgk(9)/0.030792300167387488891109020215229q0/
   data wgk(10)/0.034002130274329337836748795229551q0/
   data wgk(11)/0.037116271483415543560330625367620q0/
   data wgk(12)/0.040083825504032382074839284467076q0/
   data wgk(13)/0.042872845020170049476895792439495q0/
   data wgk(14)/0.045502913049921788909870584752660q0/
   data wgk(15)/0.047982537138836713906392255756915q0/
   data wgk(16)/0.050277679080715671963325259433440q0/
   data wgk(17)/0.052362885806407475864366712137873q0/
   data wgk(18)/0.054251129888545490144543370459876q0/
   data wgk(19)/0.055950811220412317308240686382747q0/
   data wgk(20)/0.057437116361567832853582693939506q0/
   data wgk(21)/0.058689680022394207961974175856788q0/
   data wgk(22)/0.059720340324174059979099291932562q0/
   data wgk(23)/0.060539455376045862945360267517565q0/
   data wgk(24)/0.061128509717053048305859030416293q0/
   data wgk(25)/0.061471189871425316661544131965264q0/
   !       note: wgk (26) was calculated from the values of wgk(1..25)
   data wgk(26)/0.061580818067832935078759824240066q0/

   epmach = epsilon(a)
   uflow = tiny(a)

   centr = 0.5q0*(a + b)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)
   fc = f(tx, ty, centr)
   resg = wg(13)*fc
   resk = wgk(26)*fc
   resabs = abs(resk)
   do j = 1, 12
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 13
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(26)*abs(fc - reskh)
   do j = 1, 25
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)

   return
end subroutine dqk51_3f

subroutine dqk61_3f(f, a, b, result, abserr, resabs, resasc, ty, tx)
   implicit none
   real*16 a, absc, abserr, b, centr, abs, dhlgth, dmax1, dmin1, &
      epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, resasc, &
      resg, resk, reskh, result, uflow, wg, wgk, xgk, tx, ty
   integer j, jtw, jtwm1
   external f
   dimension fv1(30), fv2(30), xgk(31), wgk(31), wg(15)

   data wg(1)/0.007968192496166605615465883474674q0/
   data wg(2)/0.018466468311090959142302131912047q0/
   data wg(3)/0.028784707883323369349719179611292q0/
   data wg(4)/0.038799192569627049596801936446348q0/
   data wg(5)/0.048402672830594052902938140422808q0/
   data wg(6)/0.057493156217619066481721689402056q0/
   data wg(7)/0.065974229882180495128128515115962q0/
   data wg(8)/0.073755974737705206268243850022191q0/
   data wg(9)/0.080755895229420215354694938460530q0/
   data wg(10)/0.086899787201082979802387530715126q0/
   data wg(11)/0.092122522237786128717632707087619q0/
   data wg(12)/0.096368737174644259639468626351810q0/
   data wg(13)/0.099593420586795267062780282103569q0/
   data wg(14)/0.101762389748405504596428952168554q0/
   data wg(15)/0.102852652893558840341285636705415q0/
   !
   data xgk(1)/0.999484410050490637571325895705811q0/
   data xgk(2)/0.996893484074649540271630050918695q0/
   data xgk(3)/0.991630996870404594858628366109486q0/
   data xgk(4)/0.983668123279747209970032581605663q0/
   data xgk(5)/0.973116322501126268374693868423707q0/
   data xgk(6)/0.960021864968307512216871025581798q0/
   data xgk(7)/0.944374444748559979415831324037439q0/
   data xgk(8)/0.926200047429274325879324277080474q0/
   data xgk(9)/0.905573307699907798546522558925958q0/
   data xgk(10)/0.882560535792052681543116462530226q0/
   data xgk(11)/0.857205233546061098958658510658944q0/
   data xgk(12)/0.829565762382768397442898119732502q0/
   data xgk(13)/0.799727835821839083013668942322683q0/
   data xgk(14)/0.767777432104826194917977340974503q0/
   data xgk(15)/0.733790062453226804726171131369528q0/
   data xgk(16)/0.697850494793315796932292388026640q0/
   data xgk(17)/0.660061064126626961370053668149271q0/
   data xgk(18)/0.620526182989242861140477556431189q0/
   data xgk(19)/0.579345235826361691756024932172540q0/
   data xgk(20)/0.536624148142019899264169793311073q0/
   data xgk(21)/0.492480467861778574993693061207709q0/
   data xgk(22)/0.447033769538089176780609900322854q0/
   data xgk(23)/0.400401254830394392535476211542661q0/
   data xgk(24)/0.352704725530878113471037207089374q0/
   data xgk(25)/0.304073202273625077372677107199257q0/
   data xgk(26)/0.254636926167889846439805129817805q0/
   data xgk(27)/0.204525116682309891438957671002025q0/
   data xgk(28)/0.153869913608583546963794672743256q0/
   data xgk(29)/0.102806937966737030147096751318001q0/
   data xgk(30)/0.051471842555317695833025213166723q0/
   data xgk(31)/0.000000000000000000000000000000000q0/
   !
   data wgk(1)/0.001389013698677007624551591226760q0/
   data wgk(2)/0.003890461127099884051267201844516q0/
   data wgk(3)/0.006630703915931292173319826369750q0/
   data wgk(4)/0.009273279659517763428441146892024q0/
   data wgk(5)/0.011823015253496341742232898853251q0/
   data wgk(6)/0.014369729507045804812451432443580q0/
   data wgk(7)/0.016920889189053272627572289420322q0/
   data wgk(8)/0.019414141193942381173408951050128q0/
   data wgk(9)/0.021828035821609192297167485738339q0/
   data wgk(10)/0.024191162078080601365686370725232q0/
   data wgk(11)/0.026509954882333101610601709335075q0/
   data wgk(12)/0.028754048765041292843978785354334q0/
   data wgk(13)/0.030907257562387762472884252943092q0/
   data wgk(14)/0.032981447057483726031814191016854q0/
   data wgk(15)/0.034979338028060024137499670731468q0/
   data wgk(16)/0.036882364651821229223911065617136q0/
   data wgk(17)/0.038678945624727592950348651532281q0/
   data wgk(18)/0.040374538951535959111995279752468q0/
   data wgk(19)/0.041969810215164246147147541285970q0/
   data wgk(20)/0.043452539701356069316831728117073q0/
   data wgk(21)/0.044814800133162663192355551616723q0/
   data wgk(22)/0.046059238271006988116271735559374q0/
   data wgk(23)/0.047185546569299153945261478181099q0/
   data wgk(24)/0.048185861757087129140779492298305q0/
   data wgk(25)/0.049055434555029778887528165367238q0/
   data wgk(26)/0.049795683427074206357811569379942q0/
   data wgk(27)/0.050405921402782346840893085653585q0/
   data wgk(28)/0.050881795898749606492297473049805q0/
   data wgk(29)/0.051221547849258772170656282604944q0/
   data wgk(30)/0.051426128537459025933862879215781q0/
   data wgk(31)/0.051494729429451567558340433647099q0/

   epmach = epsilon(a)
   uflow = tiny(a)
   !
   centr = 0.5q0*(b + a)
   hlgth = 0.5q0*(b - a)
   dhlgth = abs(hlgth)

   resg = 0.0q0
   fc = f(tx, ty, centr)
   resk = wgk(31)*fc
   resabs = abs(resk)
   do j = 1, 15
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1 + fval2
      resg = resg + wg(j)*fsum
      resk = resk + wgk(jtw)*fsum
      resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
   end do
   do j = 1, 15
      jtwm1 = j*2 - 1
      absc = hlgth*xgk(jtwm1)
      fval1 = f(tx, ty, centr - absc)
      fval2 = f(tx, ty, centr + absc)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1 + fval2
      resk = resk + wgk(jtwm1)*fsum
      resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
   end do
   reskh = resk*0.5q0
   resasc = wgk(31)*abs(fc - reskh)
   do j = 1, 30
      resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
   end do
   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk - resg)*hlgth)
   if (resasc .ne. 0.0q0 .and. abserr .ne. 0.0q0) &
      abserr = resasc*dmin1(0.1d+01, (0.2d+03*abserr/resasc)**1.5q0)
   if (resabs .gt. uflow/(0.5d+02*epmach)) abserr = dmax1 &
                                                    ((epmach*0.5d+02)*resabs, abserr)
   return
end subroutine dqk61_3f

!######

subroutine dqpsrt(limit, last, maxerr, ermax, elist, iord, nrmax)
   implicit none
   real*16 elist, ermax, errmax, errmin
   integer i, ibeg, ido, iord, isucc, j, jbnd, jupbn, k, last, limit, maxerr, &
      nrmax, igt, igt1, igt2
   dimension elist(last), iord(last)

   igt = 0
   if (last .le. 2) then
      iord(1) = 1
      iord(2) = 2
      igt = 90
   end if
   if (igt .ne. 90) then
      errmax = elist(maxerr)
      if (nrmax .ne. 1) then
         ido = nrmax - 1
         do i = 1, ido
            isucc = iord(nrmax - 1)
            if (errmax .le. elist(isucc)) exit
            iord(nrmax) = isucc
            nrmax = nrmax - 1
         end do
      end if
      jupbn = last
      if (last .gt. (limit/2 + 2)) jupbn = limit + 3 - last
      errmin = elist(last)
      jbnd = jupbn - 1
      ibeg = nrmax + 1
      igt1 = 0
      if (ibeg .le. jbnd) then
         do i = ibeg, jbnd
            isucc = iord(i)
            if (errmax .ge. elist(isucc)) then
               igt1 = 60
               exit
            end if
            iord(i - 1) = isucc
         end do
      end if
      if (igt1 .ne. 60) then
         iord(jbnd) = maxerr
         iord(jupbn) = last
         igt = 90
      end if
      igt2 = 0
      if (igt .ne. 90) then
         iord(i - 1) = maxerr
         k = jbnd
         do j = i, jbnd
            isucc = iord(k)
            if (errmin .lt. elist(isucc)) then
               igt2 = 80
               exit
            end if
            iord(k + 1) = isucc
            k = k - 1
         end do
         if (igt2 .ne. 80) then
            iord(i) = last
            igt = 90
         end if
         if (igt .ne. 90) then
            iord(k + 1) = last
         end if
      end if
   end if
   maxerr = iord(nrmax)
   ermax = elist(maxerr)
   return
end subroutine dqpsrt

subroutine xerror(xmess, nmess, nerr, level)
   !*****************************************************************************80
   !
  !! XERROR replaces the SLATEC XERROR routine.
   !
   !  Modified:
   !
   !    12 September 2015
   !
   implicit none

   integer::level, nerr, nmess
   character(len=*) xmess

   if (1 <= LEVEL) then
      WRITE (*, '(1X,A)') XMESS(1:NMESS)
      WRITE (*, '('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)') &
         NERR, LEVEL
   end if

   return
end subroutine xerror
