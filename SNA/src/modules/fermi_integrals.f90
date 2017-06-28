!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SRO_EOS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE Fermi_Integrals_Mod

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

CONTAINS
!#############################################################################
!
!   FUNCTIONs to calculate Fermi integrals of orders -1/2, 1/2, 3/2 and 5/2
!    use prescription of Toshio Fukushima
!
!    Applied Mathematics and Computation
!    Volume 259, 15 May 2015, Pages 708-
!
!   FUNCTION to calculate inverse Fermi integrals of order 1/2
!    uses the prescription of Toshio Fukushima
!
!    Applied Mathematics and Computation
!    Volume 259, 15 May 2015, Pages 698-707
!
!   FUNCTIONs are called:
!
!   F_-1/2(x) = fermi_minus_one_half(x)
!   F_+1/2(x) = fermi_one_half(x)
!   F_+3/2(x) = fermi_three_halves(x)
!   F_+5/2(x) = fermi_five_halves(x)
!
!   X_-1/2(x) = inverse_fermi_minus_one_half(x)
!   X_+1/2(x) = inverse_fermi_one_half(x)
!   X_+3/2(x) = inverse_fermi_three_halves(x)
!   X_+5/2(x) = inverse_fermi_five_halves(x)
!
!#############################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP) FUNCTION fermi_minus_three_halves(x)
!
! double precision rational minimax approximation of
!  Fermi-Dirac integral of order k=-3/2
!
! Reference: Fukushima, T. (2014, submitted to App. Math. Comp.)
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!  slightly modified by Andre da Silva Schneider 10/26/2015
!   to fit SRO_EOS source code
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: ex,t,w,s,fd
    REAL(DP), PARAMETER :: factor=-2.d0    ! = 1/(k+1)
    !
    !write(*,"(a20,1pe15.7)") "(fdm3h) x=",x
    IF(x.lt.-2.d0) THEN
        ex=exp(x)
        t=ex*7.38905609893065023d0
        fd=ex*(-3.54490770181103205d0 &
        +ex*(82737.595643818605d0 &
        +t*(18481.5553495836940d0 &
        +t*(1272.73919064487495d0 &
        +t*(26.3420403338352574d0 &
        -t*0.00110648970639283347d0 &
        ))))/(16503.7625405383183d0 &
        +t*(6422.0552658801394d0 &
        +t*(890.85389683932154d0 &
        +t*(51.251447078851450d0 &
        +t)))))
    ELSEIF(x.lt.0.d0) THEN
        s=-0.5d0*x
        t=1.d0-s
        fd=-(946.638483706348559d0 &
        +t*(76.3328330396778450d0 &
        +t*(62.7809183134124193d0 &
        +t*(83.8442376534073219d0 &
        +t*(23.2285755924515097d0 &
        +t*(3.21516808559640925d0 &
        +t*(1.58754232369392539d0 &
        +t*(0.687397326417193593d0 &
        +t*0.111510355441975495d0 &
        ))))))))/(889.4123665319664d0 &
        +s*(126.7054690302768d0 &
        +s*(881.4713137175090d0 &
        +s*(108.2557767973694d0 &
        +s*(289.38131234794585d0 &
        +s*(27.75902071820822d0 &
        +s*(34.252606975067480d0 &
        +s*(1.9592981990370705d0 &
        +s))))))))
    ELSEIF(x.lt.2.d0) THEN
        t=0.5d0*x
        fd=-(754.61690882095729d0 &
        +t*(565.56180911009650d0 &
        +t*(494.901267018948095d0 &
        +t*(267.922900418996927d0 &
        +t*(110.418683240337860d0 &
        +t*(39.4050164908951420d0 &
        +t*(10.8654460206463482d0 &
        +t*(2.11194887477009033d0 &
        +t*0.246843599687496060d0 &
        ))))))))/(560.03894899770103d0 &
        +t*(70.007586553114572d0 &
        +t*(582.42052644718871d0 &
        +t*(56.181678606544951d0 &
        +t*(205.248662395572799d0 &
        +t*(12.5169006932790528d0 &
        +t*(27.2916998671096202d0 &
        +t*(0.53299717876883183d0 &
        +t))))))))
    ELSEIF(x.lt.5.d0) THEN
        t=0.3333333333333333333d0*(x-2.d0)
        fd=-(526.022770226139287d0 &
        +t*(631.116211478274904d0 &
        +t*(516.367876532501329d0 &
        +t*(267.894697896892166d0 &
        +t*(91.3331816844847913d0 &
        +t*(17.5723541971644845d0 &
        +t*(1.46434478819185576d0 &
        +t*(1.29615441010250662d0 &
        +t*0.223495452221465265d0 &
        ))))))))/(354.867400305615304d0 &
        +t*(560.931137013002977d0 &
        +t*(666.070260050472570d0 &
        +t*(363.745894096653220d0 &
        +t*(172.272943258816724d0 &
        +t*(23.7751062504377332d0 &
        +t*(12.5916012142616255d0 &
        +t*(-0.888604976123420661d0 &
        +t))))))))
    ELSEIF(x.lt.10.d0) THEN
        t=0.2d0*x-1.d0
        fd=-(18.0110784494455205d0 &
        +t*(36.1225408181257913d0 &
        +t*(38.4464752521373310d0 &
        +t*(24.1477896166966673d0 &
        +t*(9.27772356782901602d0 &
        +t*(2.49074754470533706d0 &
        +t*(0.163824586249464178d0 &
        -t*0.00329391807590771789d0 &
        )))))))/(18.8976860386360201d0 &
        +t*(49.3696375710309920d0 &
        +t*(60.9273314194720251d0 &
        +t*(43.6334649971575003d0 &
        +t*(20.6568810936423065d0 &
        +t*(6.11094689399482273d0 &
        +t))))))
    ELSEIF(x.lt.20.d0) THEN
        t=0.1d0*x-1.d0
        fd=-(4.10698092142661427d0 &
        +t*(17.1412152818912658d0 &
        +t*(32.6347877674122945d0 &
        +t*(36.6653101837618939d0 &
        +t*(25.9424894559624544d0 &
        +t*(11.2179995003884922d0 &
        +t*(2.30099511642112478d0 &
        +t*(0.0928307248942099967d0 &
        -t*0.00146397877054988411d0 &
        ))))))))/(6.40341731836622598d0 &
        +t*(30.1333068545276116d0 &
        +t*(64.0494725642004179d0 &
        +t*(80.5635003792282196d0 &
        +t*(64.9297873014508805d0 &
        +t*(33.3013900893183129d0 &
        +t*(9.61549304470339929d0 &
        +t)))))))
    ELSEIF(x.lt.40.d0) THEN
        t=0.05d0*x-1.d0
        fd=-(95.2141371910496454d0 &
        +t*(420.050572604265456d0 &
        +t*(797.778374374075796d0 &
        +t*(750.378359146985564d0 &
        +t*(324.818150247463736d0 &
        +t*(50.3115388695905757d0 &
        +t*(0.372431961605507103d0 &
        +t*(-0.103162211894757911d0 &
        +t*0.00191752611445211151d0 &
        ))))))))/(212.232981736099697d0 &
        +t*(1043.79079070035083d0 &
        +t*(2224.50099218470684d0 &
        +t*(2464.84669868672670d0 &
        +t*(1392.55318009810070d0 &
        +t*(346.597189642259199d0 &
        +t*(22.7314613168652593d0 &
        -t)))))))
    ELSE
        w=1.d0/(x*x)
        s=1.d0-1600.d0*w
        fd=factor/sqrt(x)*(1.d0 &
        +w*(12264.3569103180524d0 &
        +s*(3204.34872454052352d0 &
        +s*(140.119604748253961d0 &
        +s*0.523918919699235590d0 &
        )))/(9877.87829948067200d0 &
        +s*(2644.71979353906092d0 &
        +s*(128.863768007644572d0 &
        +s))))
    ENDIF

    fermi_minus_three_halves=fd

END FUNCTION fermi_minus_three_halves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  REAL(DP) FUNCTION fermi_minus_one_half(x)
!
! double precision rational minimax approximation
!   of Fermi-Dirac integral of order k=-1/2
!
! Reference: Fukushima, T. (2014, submitted to App. Math. Comp.)
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!  slightly modIFied by Andre da Silva Schneider 10/26/2015
!   to fit SRO_EOS source code
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: ex,t,w,s,fd
    REAL(DP), PARAMETER :: factor=2.d0    ! = 1/(k+1)

    IF(x.lt.-2.d0) THEN
      ex=exp(x)
      t=ex*7.38905609893065023d0
      fd=ex*(1.77245385090551603d0 &
      -ex*(40641.4537510284430d0 &
      +t*(9395.7080940846442d0 &
      +t*(649.96168315267301d0 &
      +t*(12.7972295804758967d0 &
      +t*0.00153864350767585460d0 &
      ))))/(32427.1884765292940d0 &
      +t*(11079.9205661274782d0 &
      +t*(1322.96627001478859d0 &
      +t*(63.738361029333467d0 &
      +t)))))
    ELSEIF(x.lt.0.d0) THEN
      s=-0.5d0*x
      t=1.d0-s
      fd=(272.770092131932696d0 &
      +t*(30.8845653844682850d0 &
      +t*(-6.43537632380366113d0 &
      +t*(14.8747473098217879d0 &
      +t*(4.86928862842142635d0 &
      +t*(-1.53265834550673654d0 &
      +t*(-1.02698898315597491d0 &
      +t*(-0.177686820928605932d0 &
      -t*0.00377141325509246441d0 &
      ))))))))/(293.075378187667857d0 &
      +s*(305.818162686270816d0 &
      +s*(299.962395449297620d0 &
      +s*(207.640834087494249d0 &
      +s*(92.0384803181851755d0 &
      +s*(37.0164914112791209d0 &
      +s*(7.88500950271420583d0 &
      +s)))))))
    ELSEIF(x.lt.2.d0) THEN
      t=0.5d0*x
      fd=(3531.50360568243046d0 &
      +t*(6077.5339658420037d0 &
      +t*(6199.7700433981326d0 &
      +t*(4412.78701919567594d0 &
      +t*(2252.27343092810898d0 &
      +t*(811.84098649224085d0 &
      +t*(191.836401053637121d0 &
      +t*23.2881838959183802d0 &
      )))))))/(3293.83702584796268d0 &
      +t*(1528.97474029789098d0 &
      +t*(2568.48562814986046d0 &
      +t*(925.64264653555825d0 &
      +t*(574.23248354035988d0 &
      +t*(132.803859320667262d0 &
      +t*(29.8447166552102115d0 &
      +t)))))))
    ELSEIF(x.lt.5.d0) THEN
      t=0.3333333333333333333d0*(x-2.d0)
      fd=(4060.70753404118265d0 &
      +t*(10812.7291333052766d0 &
      +t*(13897.5649482242583d0 &
      +t*(10628.4749852740029d0 &
      +t*(5107.70670190679021d0 &
      +t*(1540.84330126003381d0 &
      +t*(284.452720112970331d0 &
      +t*29.5214417358484151d0 &
      )))))))/(1564.58195612633534d0 &
      +t*(2825.75172277850406d0 &
      +t*(3189.16066169981562d0 &
      +t*(1955.03979069032571d0 &
      +t*(828.000333691814748d0 &
      +t*(181.498111089518376d0 &
      +t*(32.0352857794803750d0 &
      +t)))))))
    ELSEIF(x.lt.10.d0) THEN
      t=0.2d0*x-1.d0
      fd=(1198.41719029557508d0 &
      +t*(3263.51454554908654d0 &
      +t*(3874.97588471376487d0 &
      +t*(2623.13060317199813d0 &
      +t*(1100.41355637121217d0 &
      +t*(267.469532490503605d0 &
      +t*(25.4207671812718340d0 &
      +t*0.389887754234555773d0 &
      )))))))/(273.407957792556998d0 &
      +t*(595.918318952058643d0 &
      +t*(605.202452261660849d0 &
      +t*(343.183302735619981d0 &
      +t*(122.187622015695729d0 &
      +t*(20.9016359079855933d0 &
      +t))))))
    ELSEIF(x.lt.20.d0) THEN
      t=0.1d0*x-1.d0
      fd=(9446.00169435237637d0 &
      +t*(36843.4448474028632d0 &
      +t*(63710.1115419926191d0 &
      +t*(62985.2197361074768d0 &
      +t*(37634.5231395700921d0 &
      +t*(12810.9898627807754d0 &
      +t*(1981.56896138920963d0 &
      +t*81.4930171897667580d0 &
      )))))))/(1500.04697810133666d0 &
      +t*(5086.91381052794059d0 &
      +t*(7730.01593747621895d0 &
      +t*(6640.83376239360596d0 &
      +t*(3338.99590300826393d0 &
      +t*(860.499043886802984d0 &
      +t*(78.8565824186926692d0 &
      +t)))))))
    ELSEIF(x.lt.40.d0) THEN
      t=0.05d0*x-1.d0
      fd=(22977.9657855367223d0 &
      +t*(123416.616813887781d0 &
      +t*(261153.765172355107d0 &
      +t*(274618.894514095795d0 &
      +t*(149710.718389924860d0 &
      +t*(40129.3371700184546d0 &
      +t*(4470.46495881415076d0 &
      +t*132.684346831002976d0 &
      )))))))/(2571.68842525335676d0 &
      +t*(12521.4982290775358d0 &
      +t*(23268.1574325055341d0 &
      +t*(20477.2320119758141d0 &
      +t*(8726.52577962268114d0 &
      +t*(1647.42896896769909d0 &
      +t*(106.475275142076623d0 &
      +t)))))))
    ELSE
      w=1.d0/(x*x)
      t=1600.d0*w
      fd=sqrt(x)*factor*(1.d0 &
      -w*(0.411233516712009968d0 &
      +t*(0.00110980410034088951d0 &
      +t*(0.0000113689298990173683d0 &
      +t*(2.56931790679436797d-7 &
      +t*(9.97897786755446178d-9 &
      +t*8.67667698791108582d-10))))))
    ENDIF

    fermi_minus_one_half = fd

  END FUNCTION fermi_minus_one_half
!#############################################################################

  REAL(DP) FUNCTION fermi_one_half(x)

    IMPLICIT NONE
!
! double precision rational minimax approximation
!  of Fermi-Dirac integral of order k=1/2
!
! Reference: Fukushima, T. (2014, submitted to App. Math. Comp.)
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!  slightly modIFied by Andre da Silva Schneider 10/26/2015
!   to fit LSEOSF90 source code

    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: ex,t,w,s,fd
    REAL(DP), PARAMETER :: factor=2.d0/3.d0    ! = 1/(k+1)

    IF(x.lt.-2.d0) THEN
      ex=exp(x)
      t=ex*7.38905609893065023d0
      fd=ex*(0.886226925452758014d0 &
      -ex*(19894.4553386951666d0 &
      +t*(4509.64329955948557d0 &
      +t*(303.461789035142376d0 &
      +t*(5.7574879114754736d0 &
      +t*0.00275088986849762610d0 &
      ))))/(63493.915041308052d0 &
      +t*(19070.1178243603945d0 &
      +t*(1962.19362141235102d0 &
      +t*(79.250704958640158d0 &
      +t)))))
    ELSEIF(x.lt.0.d0) THEN
      s=-0.5d0*x
      t=1.d0-s
      fd=(149.462587768865243d0 &
      +t*(22.8125889885050154d0 &
      +t*(-0.629256395534285422d0 &
      +t*(9.08120441515995244d0 &
      +t*(3.35357478401835299d0 &
      +t*(-0.473677696915555805d0 &
      +t*(-0.467190913556185953d0 &
      +t*(-0.0880610317272330793d0 &
      -t*0.00262208080491572673d0 &
      ))))))))/(269.94660938022644d0 &
      +s*(343.6419926336247d0 &
      +s*(323.9049470901941d0 &
      +s*(218.89170769294024d0 &
      +s*(102.31331350098315d0 &
      +s*(36.319337289702664d0 &
      +s*(8.3317401231389461d0 &
      +s)))))))
    ELSEIF(x.lt.2.d0) THEN
      t=0.5d0*x
      fd=(71652.717119215557d0 &
      +t*(134954.734070223743d0 &
      +t*(153693.833350315645d0 &
      +t*(123247.280745703400d0 &
      +t*(72886.293647930726d0 &
      +t*(32081.2499422362952d0 &
      +t*(10210.9967337762918d0 &
      +t*(2152.71110381320778d0 &
      +t*232.906588165205042d0 &
      ))))))))/(105667.839854298798d0 &
      +t*(31946.0752989314444d0 &
      +t*(71158.788776422211d0 &
      +t*(15650.8990138187414d0 &
      +t*(13521.8033657783433d0 &
      +t*(1646.98258283527892d0 &
      +t*(618.90691969249409d0 &
      +t*(-3.36319591755394735d0 &
      +t))))))))
    ELSEIF(x.lt.5.d0) THEN
      t=0.3333333333333333333d0*(x-2.d0)
      fd=(23744.8706993314289d0 &
      +t*(68257.8589855623002d0 &
      +t*(89327.4467683334597d0 &
      +t*(62766.3415600442563d0 &
      +t*(20093.6622609901994d0 &
      +t*(-2213.89084119777949d0 &
      +t*(-3901.66057267577389d0 &
      -t*948.642895944858861d0 &
      )))))))/(9488.61972919565851d0 &
      +t*(12514.8125526953073d0 &
      +t*(9903.44088207450946d0 &
      +t*(2138.15420910334305d0 &
      +t*(-528.394863730838233d0 &
      +t*(-661.033633995449691d0 &
      +t*(-51.4481470250962337d0 &
      +t)))))))
    ELSEIF(x.lt.10.d0) THEN
      t=0.2d0*x-1.d0
      fd=(311337.452661582536d0 &
      +t*(1.11267074416648198d6 &
      +t*(1.75638628895671735d6 &
      +t*(1.59630855803772449d6 &
      +t*(910818.935456183774d0 &
      +t*(326492.733550701245d0 &
      +t*(65507.2624972852908d0 &
      +t*4809.45649527286889d0 &
      )))))))/(39721.6641625089685d0 &
      +t*(86424.7529107662431d0 &
      +t*(88163.7255252151780d0 &
      +t*(50615.7363511157353d0 &
      +t*(17334.9774805008209d0 &
      +t*(2712.13170809042550d0 &
      +t*(82.2205828354629102d0 &
      -t)))))))*0.999999999999999877d0
    ELSEIF(x.lt.20.d0) THEN
      t=0.1d0*x-1.d0
      fd=(7.26870063003059784d6 &
      +t*(2.79049734854776025d7 &
      +t*(4.42791767759742390d7 &
      +t*(3.63735017512363365d7 &
      +t*(1.55766342463679795d7 &
      +t*(2.97469357085299505d6 &
      +t*154516.447031598403d0 &
      ))))))/(340542.544360209743d0 &
      +t*(805021.468647620047d0 &
      +t*(759088.235455002605d0 &
      +t*(304686.671371640343d0 &
      +t*(39289.4061400542309d0 &
      +t*(582.426138126398363d0 &
      +t*(11.2728194581586028d0 &
      -t)))))))
    ELSEIF(x.lt.40.d0) THEN
      t=0.05d0*x-1.d0
      fd=(4.81449797541963104d6 &
      +t*(1.85162850713127602d7 &
      +t*(2.77630967522574435d7 &
      +t*(2.03275937688070624d7 &
      +t*(7.41578871589369361d6 &
      +t*(1.21193113596189034d6 &
      +t*63211.9545144644852d0 &
      ))))))/(80492.7765975237449d0 &
      +t*(189328.678152654840d0 &
      +t*(151155.890651482570d0 &
      +t*(48146.3242253837259d0 &
      +t*(5407.08878394180588d0 &
      +t*(112.195044410775577d0 &
      -t))))))
    ELSE
      w=1.d0/(x*x)
      s=1.d0-1600.d0*w
      fd=x*sqrt(x)*0.666666666666666667d0*(1.d0+w &
      *(8109.79390744477921d0 &
      +s*(342.069867454704106d0 &
      +s*1.07141702293504595d0)) &
      /(6569.98472532829094d0 &
      +s*(280.706465851683809d0 &
      +s)))
    ENDIF

    fermi_one_half = fd

  END FUNCTION fermi_one_half

!#############################################################################

  REAL(DP) FUNCTION fermi_three_halves(x)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: ex,t,w,s,fd
    REAL(DP), PARAMETER :: factor=2.d0/5.d0    ! = 1/(k+1)

! double precision rational minimax approximation of Fermi-Dirac integral of order k=3/2
!
! Reference: Fukushima, T. (2014, submitted to App. Math. Comp.)
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
!
!write(*,"(a20,1pe15.7)") "(fd3h) x=",x
    IF(x.lt.-2.d0) THEN
      ex=exp(x)
      t=ex*7.38905609893065023d0
      fd=ex*(1.32934038817913702d0 &
      -ex*(1346.14119046566636d0 &
      +t*(199.946876779712712d0 &
      +t*(6.5210149677288048d0 &
      +t*0.0108588591982722183d0 &
      )))/(5728.3481201778541d0 &
      +t*(1132.17837281710987d0 &
      +t*(64.805243148002602d0 &
      +t))))
    ELSEIF(x.lt.0.d0) THEN
      s=-0.5d0*x
      t=1.d0-s
      fd=(631.667081787115831d0 &
      +t*(504.131655805666135d0 &
      +t*(113.449065431934917d0 &
      +t*(56.0939647772947784d0 &
      +t*(43.3374223200846752d0 &
      +t*(12.8047010597109577d0 &
      +t*(0.219164386586949410d0 &
      +t*(-0.678659552658390139d0 &
      -t*0.126533769309899232d0 &
      ))))))))/(1180.5112183558028d0 &
      +s*(1101.0159189871135d0 &
      +s*(864.4448234404281d0 &
      +s*(392.2018227840790d0 &
      +s*(89.58093202779063d0 &
      +s*(-9.95066218572899d0 &
      +s*(-17.312068771997626d0 &
      +s*(-6.4116162917822773d0 &
      -s))))))))
    ELSEIF(x.lt.2.d0) THEN
      t=0.5d0*x
      fd=(90122.488639370400d0 &
      +t*(157095.208147064037d0 &
      +t*(166879.962599668589d0 &
      +t*(125708.597728045460d0 &
      +t*(69968.278213181390d0 &
      +t*(29035.3292989055404d0 &
      +t*(8736.4439472398517d0 &
      +t*(1747.16784760309227d0 &
      +t*180.132410666734053d0 &
      ))))))))/(78176.777123671727d0 &
      +t*(-1681.44633240543085d0 &
      +t*(38665.7913035496031d0 &
      +t*(-2527.29685826087874d0 &
      +t*(5062.6683078100048d0 &
      +t*(-553.21165462054589d0 &
      +t*(165.395637981775430d0 &
      +t*(- 18.0295465153725544d0 &
      +t))))))))
    ELSEIF(x.lt.5.d0) THEN
      t=0.3333333333333333333d0*(x-2.d0)
      fd=(912944.432058014054d0 &
      +t*(3.28217091334054338d6 &
      +t*(5.59250227196369585d6 &
      +t*(5.76136129685687470d6 &
      +t*(3.84331519034749983d6 &
      +t*(1.65284168824947710d6 &
      +t*(423452.676670436605d0 &
      +t*49835.4127241373113d0 &
      )))))))/(164873.145721762182d0 &
      +t*(257442.511191094986d0 &
      +t*(225604.160532840884d0 &
      +t*(99932.1955662320024d0 &
      +t*(24761.0878784286761d0 &
      +t*(1398.26392212830777d0 &
      +t*(-36.4450237523474167d0 &
      +t)))))))
    ELSEIF(x.lt.10.d0) THEN
      t=0.2d0*x-1.d0
      fd=(1.88412548327216052d6 &
      +t*(8.08838896259910792d6 &
      +t*(1.56869793001790529d7 &
      +t*(1.79109792599373447d7 &
      +t*(1.31345142328147214d7 &
      +t*(6.29500412046744325d6 &
      +t*(1.89326213154091054d6 &
      +t*(312372.643127575407d0 &
      +t*18814.7420442630170d0 &
      ))))))))/(67768.3347951202583d0 &
      +t*(147635.914444221358d0 &
      +t*(151908.303165069423d0 &
      +t*(86671.1222110642970d0 &
      +t*(27855.9481608626219d0 &
      +t*(3833.22697473114940d0 &
      +t*(98.3384567064269554d0 &
      -t)))))))*0.999999999999999876d0  ! correction to remove bias
    ELSEIF(x.lt.20.d0) THEN
      t=0.1d0*x-1.d0
      fd=(1.59656593348660977d9 &
      +t*(7.32769737561517060d9 &
      +t*(1.42662658588280191d10 &
      +t*(1.51238422045169918d10 &
      +t*(9.27233604548095476d9 &
      +t*(3.18834513406577423d9 &
      +t*(5.36061988605886123d8 &
      +t*3.03619219668246382d7 &
      )))))))/(1.18906980815759995d7 &
      +t*(2.62209219322122975d7 &
      +t*(2.28143701746618729d7 &
      +t*(8.57156701742181874d6 &
      +t*(1.13860063870524239d6 &
      +t*(27091.7884208687379d0 &
      +t*(-275.664733379090447d0 &
      +t)))))))*0.999999999999999829d0  ! correction to remove bias
    ELSEIF(x.lt.40.d0) THEN
      t=0.05d0*x-1.d0
      fd=(2.60437581212904589d8 &
      +t*(1.08771546307370080d9 &
      +t*(1.81531350939088943d9 &
      +t*(1.52833764636304939d9 &
      +t*(6.70684451492750149d8 &
      +t*(1.40870639531414149d8 &
      +t*1.04957900377463854d7 &
      ))))))/(358448.871166784200d0 &
      +t*(611808.419702466190d0 &
      +t*(326307.561591723775d0 &
      +t*(58407.9904827573816d0 &
      +t*(2049.50040323021794d0 &
      +t*(-39.8767861209088081d0 &
      +t))))))*0.999999999999999828d0
    ELSE
      w=1.d0/(x*x)
      s=1.d0-1600.d0*w
      fd=x*x*sqrt(x)*factor*(1.d0 &
      +w*(6.16739021212286242d0 &
      +s*(0.00111530123694574981d0 &
      +s*(-2.79156524536560815d-6 &
      +s*(2.95571462110856359d-8 &
      -s*6.70917556862133933d-10)))))
    ENDIF

    fermi_three_halves=fd

  END FUNCTION fermi_three_halves


  REAL(DP) FUNCTION inverse_fermi_one_half(x)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: s, t, u, v, w, y, z
!
! double precision rational minimax approximation of
! inverse Fermi-Dirac integral of order k=1/2
!
! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 698-707)
!
! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!
    u = x

    IF(u.lt.1.17683303804380831d0) THEN
      t=u*0.849738210666018375d0
      z=t*(156377.8333056294d0 &
      +t*(48177.5705898287d0 &
      +t*(5847.07218383812d0 &
      +t*(335.3978079672194d0 &
      +t*7.84411868029912d0 &
      ))))/(117762.02905535089d0 &
      +t*(-19007.26938370368d0 &
      +t*(1376.2936928453140d0 &
      +t*(-54.11372698481717d0 &
      +t))))
      y=log(z)
    ELSEIF(u.lt.3.82993088157949761d0) THEN
      t=0.376917874490198033d0*u &
      -0.443569407329314587d0
      y=(489.140447310410217d0 &
      +t*(5335.07269317261966d0 &
      +t*(20169.0736140442509d0 &
      +t*(35247.8115595510907d0 &
      +t*(30462.3668614714761d0 &
      +t*(12567.9032426128967d0 &
      +t*(2131.86789357398657d0 &
      +t*93.6520172085419439d0 &
      )))))))/(656.826207643060606d0 &
      +t*(4274.82831051941605d0 &
      +t*(10555.7581310151498d0 &
      +t*(12341.8742094611883d0 &
      +t*(6949.18854413197094d0 &
      +t*(1692.19650634194002d0 &
      +t*(129.221772991589751d0 &
      +t)))))))
    ELSEIF(u.lt.13.3854493161866553d0) THEN
      t=0.104651569335924949d0*u &
      -0.400808277205416960d0
      y=(1019.84886406642351d0 &
      +t*(9440.18255003922075d0 &
      +t*(33947.6616363762463d0 &
      +t*(60256.7280980542786d0 &
      +t*(55243.0045063055787d0 &
      +t*(24769.8354802210838d0 &
      +t*(4511.77288617668292d0 &
      +t*211.432806336150141d0 &
      )))))))/(350.502070353586442d0 &
      +t*(2531.06296201234050d0 &
      +t*(6939.09850659439245d0 &
      +t*(9005.40197972396592d0 &
      +t*(5606.73612994134056d0 &
      +t*(1488.76634564005075d0 &
      +t*(121.537028889412581d0 &
      +t)))))))
    ELSEIF(u.lt.53.2408277860982205d0) THEN
      t=0.0250907164450825724d0*u &
      -0.335850513282463787d0
      y=(11885.8779398399498d0 &
      +t*(113220.250825178799d0 &
      +t*(408524.373881197840d0 &
      +t*(695674.357483475952d0 &
      +t*(569389.917088505552d0 &
      +t*(206433.082013681440d0 &
      +t*(27307.2535671974100d0 &
      +t*824.430826794730740d0 &
      )))))))/(1634.40491220861182d0 &
      +t*(12218.1158551884025d0 &
      +t*(32911.7869957793233d0 &
      +t*(38934.6963039399331d0 &
      +t*(20038.8358438225823d0 &
      +t*(3949.48380897796954d0 &
      +t*(215.607404890995706d0 &
      +t)))))))
    ELSEIF(u.lt.188.411871723022843d0) THEN
      t=0.00739803415638806339d0*u &
      -0.393877462475929313d0
      y=(11730.7011190435638d0 &
      +t*(99421.7455796633651d0 &
      +t*(327706.968910706902d0 &
      +t*(530425.668016563224d0 &
      +t*(438631.900516555072d0 &
      +t*(175322.855662315845d0 &
      +t*(28701.9605988813884d0 &
      +t*1258.20914464286403d0 &
      )))))))/(634.080470383026173d0 &
      +t*(4295.63159860265838d0 &
      +t*(10868.5260668911946d0 &
      +t*(12781.6871997977069d0 &
      +t*(7093.80732100760563d0 &
      +t*(1675.06417056300026d0 &
      +t*(125.750901817759662d0 &
      +t)))))))
    ELSE
      v=u**(-4.d0/3.d0)
      s=1080.13412050984017d0*v
      t=1.d0-s
      w=(1.12813495144821933d7 &
      +t*(420368.911157160874d0 &
      +t*(1689.69475714536117d0 &
      +t)))/(s*(6088.08350831295857d0 &
      +t*(221.445236759466761d0 &
      +t*0.718216708695397737d0 &
      )))
      y=sqrt(w)
    ENDIF

    inverse_fermi_one_half = y

  END FUNCTION inverse_fermi_one_half

  REAL(DP) FUNCTION fhalf(x)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x

    fhalf = fermi_one_half(x)/fermi_minus_one_half(x)

  END FUNCTION fhalf

END MODULE Fermi_Integrals_Mod
