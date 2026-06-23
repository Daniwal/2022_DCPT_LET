*======================================================================
* fluka_custom_sobp_letmom.f
*
* Combined custom FLUKA routines for the FLUKA vs SHIELD-HIT12A
* proton therapy benchmark.
*
* Contains:
*   1) SOBP source sampler
*   2) FLUSCW LET-moment scoring weights
*
* Built from canonical files:
*   custom_f/source_sampler.f
*   custom_f/fluscw_let.f
*======================================================================

*$ CREATE SOURCE.FOR
*COPY SOURCE

!> @brief
!! Particle source for pencil beam scanning in hadrontherapy.
!! It mimics the beam coming out from gantry nozzle and aims
!! at creating spread-out Bragg peak shape in depth
!! Simulated source consists of several beamlets,
!! each of them characterized by relative weight,
!! energy and size (both defined at source location).
!! Additionally energy spread may be specified by the user.
!! User can also decide whether divergent beam is emitted from
!! a reference/source plane or point-like virtual source is used.
!! This source needs an additional file (typically sobp.dat)
!! with description of pencil beam geometry and kinematics.
!!
!! This implementation is based on a template from $FLUPRO/usermvax/source.f
!!
!! In order to use the source, first compile this file using
!! following command or Flair GUI:
!!
!!  ldpmqmd -oflukadpm_sobp source_sampler.f
!!
!! Then get a file called sobp.dat and put it in the same directory as
!! your Fluka input file. In the input file add a card called SOURCE
!! to activate this custom source. To run it, call (or use Flair):
!!
!! rfluka -N0 -M1 -e flukadpm_sobp your_input_file
!!
!!
!!
!! @details
!!
!! -------------------- SOBP CONFIG FILE -------------------------------------
!! Input file (typically sobp.dat) is a text file.
!! Comment lines start with #. It may have 5,6,7,9 or 11 columns.
!! These columns can be:
!!
!!  - 5 columns: E, X, Y, FWHM, W
!!  - 6 columns: E, X, Y, FWHM_X, FWHM_Y, W
!!  - 7 columns: E, DE, X, Y, FWHM_X, FWHM_Y, W
!!  - 9 columns: E, DE, X, Y, FWHM_X, FWHM_Y, DIVX, DIVY, W
!!  - 11 columns: E, DE, X, Y, FWHM_X, FWHM_Y, DIVX, DIVY, CORX, CORY, W
!!
!! where:
!!
!!  - E : kinetic energy in GeV/nucleon
!!  - DE : energy spread (sigma) in GeV/nucleon
!!  - X, Y : position (in cm) of the beamlet/spot center
!!  - FWHM_X, FWHM_Y, FWHM : spot size (in cm)
!!  - DIVX, DIVY : angular divergence (in mrad)
!!  - CORX, CORY : correlation coefficient rho (dimensionless)
!!
!! All above-mentioned quantities are defined at beam source location (typically nozzle exit)
!!
!! For more details, see SHIELD-HIT12A manual http://shieldhit.org/index.php?id=documentation
!!
!!
!!
!! -------------------- Beam configuration in the INPUT FILE -------------------------------------
!! Beam momentum/energy, its spread and divergence is usually specified in BEAM input card.
!! Also beam width is specified there.
!! All these values except momentum spread will be ignored and overridden.
!! In case 5- or 6-columns sobp.dat file format (the case where DE is missing) is used,
!! then beam momentum will be taken from the input card, otherwise is will be ignored and overriden.
!!
!! Center of the beam spot is usually defined in BEAMPOS card in the input file,
!! which contains X,Y,Z positions and direction cosines.
!! Whatever is specified as X,Y and direction cosines will be ignored and overridden
!! by this custom source. SDUM value in BEAMPOS card will also be ignored.
!!
!! In order to activate this custom source, please add SOURCE card to the input file.
!! WHAT(1), if nonzero, enables virtual-source geometry using SADx and SADy
!! from WHAT(3) and WHAT(4). If WHAT(1) is zero or omitted, this mode is off.
!! The SOURCE card SDUM is the filename containing the spotlist. If SDUM is
!! omitted, this routine reads sobp.dat.
!! An example SOURCE card looks like that:
!!
!! SOURCE             1.0       0.0     205.0     205.0                  sobp.dat
!!
!!
!! For more details, see FLUKA documentation:
!!  - BEAM    card http://www.fluka.org/fluka.php?id=man_onl&sub=12
!!  - BEAMPOS card http://www.fluka.org/fluka.php?id=man_onl&sub=14
!!  - SOURCE  card http://www.fluka.org/fluka.php?id=man_onl&sub=71
!!
!! -------------------------------------------------------------------------------------------------
!! ==================================================================================================
!! ==================================================================================================

!> @brief
!! Returns number of text blocks separated with white characters (spaces)
!! It corresponds to number of columns in CSV file.
!! @param[in] STRING
!! @retval LENGTH number of columns
      INTEGER FUNCTION NCOLS( ASTRING )
      IMPLICIT NONE
      CHARACTER*(*) ASTRING   ! input - string of arbitrary length
      INTEGER I               !
      NCOLS = 0

*     make a loop and check for non-space followed by space
*     each occurence will mark end of a column
      DO I = 1, LEN (ASTRING) - 1
         IF ((ASTRING (I:I) .NE. ' ') .AND.
     & (ASTRING(I+1:I+1) .EQ. ' ')) THEN
            NCOLS = NCOLS + 1
         END IF

*     in a special case when last character isn't a space, add another column
      IF ( ASTRING( LEN(ASTRING):LEN(ASTRING) ) .NE. ' ' ) THEN
        NCOLS = NCOLS + 1
      END IF

      END DO

      RETURN
      END


!! ==================================================================================================
!! ==================================================================================================

!> @brief
!! Reads beam configuration file
!! @param[in] FILEPATH  path to the beam configuration file
!! @param[out] ENERGY  particle energy in GeV/nucleon
!! @param[out] DELTAE  particle energy spread (sigma) in GeV/nucleon
!! @param[out] XPOS beam spot center (X coordinate), in cm
!! @param[out] YPOS beam spot center (Y coordinate), in cm
!! @param[out] FWHMX beam spot size in X axis, in cm
!! @param[out] FWHMY beam spot size in Y axis, in cm
!! @param[out] DIVX beam spot angular divergence in X axis, in mrad
!! @param[out] DIVY beam spot angular divergence in Y axis, in mrad
!! @param[out] CORX correlation coefficient rho(x,tx) (dimensionless)
!! @param[out] CORY correlation coefficient rho(y,ty) (dimensionless)
!! @param[out] PART beamlet weight
!! @param[out] NCOLUMNS number of columns in the file, negative if file missing or corrupted
!! @param[out] NWEIGHT number of data rows in the file, negative if file missing or corrupted
      SUBROUTINE READSOBP ( FILEPATH,
     $   ENERGY, DELTAE, XPOS, YPOS, FWHMX, FWHMY,
     $   DIVX, DIVY, CORX, CORY, PART,
     $   NCOLUMNS, NWEIGHT )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'

      CHARACTER*(*) FILEPATH   ! path to the sobp file
      CHARACTER(8192) LINE
      DOUBLE PRECISION ENERGY(65000), DELTAE(65000)
      DOUBLE PRECISION XPOS(65000), YPOS(65000)
      DOUBLE PRECISION FWHMX(65000), FWHMY(65000)
      DOUBLE PRECISION DIVX(65000), DIVY(65000)
      DOUBLE PRECISION CORX(65000), CORY(65000)
      DOUBLE PRECISION PART(65000)
      INTEGER NCOLUMNS
      INTEGER NWEIGHT
      LOGICAL LEXISTS

      NWEIGHT = 0
      NCOLUMNS = 0

      WRITE(LUNOUT,*) 'SOBP SOURCE READING ', FILEPATH

*     warn if the spotlist file is missing, stop calculation
      INQUIRE( FILE=FILEPATH, EXIST=LEXISTS )
      IF ( .NOT. LEXISTS ) THEN
         WRITE(LUNOUT,*) 'SOBP FILE missing: ', FILEPATH
         RETURN
      END IF
*
*     open spotlist file for reading
      OPEN( 44, FILE=FILEPATH, STATUS='OLD' )
*
*     we will now probe the file to get the number of columns with numbers
*     we skip comment lines, first non-comment line will be saved to LINE
*     later we rewind the file to the beginning and read it with proper column format
      LINE(1:1) = '#'
      DO WHILE( LINE(1:1) .EQ. '#' )
         READ(44,'(A)',END=10) LINE
      END DO
      REWIND(44)
*     WRITE(LUNOUT,*) 'SOBP FIRST LINE', TRIM(LINE)
*
*     get number of columns from the first non-comment line
      NCOLUMNS = NCOLS(LINE)
      WRITE(LUNOUT,*) 'SOBP NUMBER OF COLUMNS', NCOLUMNS


      DO
*        fortran arrays start with 1, so we increase the counter as the loop starts
         NWEIGHT = NWEIGHT + 1
         IF (NWEIGHT .GT. 65000) THEN
            WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: too many beamlets'
            RETURN
         ENDIF

*        skip comment lines, first non-comment line will be saved to LINE
         LINE(1:1) = '#'
         DO WHILE( LINE(1:1) .EQ. '#' )
            READ(44,'(A)',END=10) LINE
         END DO

*        read the line, and guess the format from number of columns
         IF ( NCOLUMNS .EQ. 5 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), XPOS(NWEIGHT),
     $         YPOS(NWEIGHT), FWHMX(NWEIGHT), PART(NWEIGHT)
            FWHMY(NWEIGHT)  = FWHMX(NWEIGHT)
            DELTAE(NWEIGHT) = 0.0D0
            DIVX(NWEIGHT)    = 0.0D0
            DIVY(NWEIGHT)    = 0.0D0
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 6 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), XPOS(NWEIGHT),
     $         YPOS(NWEIGHT), FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         PART(NWEIGHT)
            DELTAE(NWEIGHT) = 0.0D0
            DIVX(NWEIGHT)    = 0.0D0
            DIVY(NWEIGHT)    = 0.0D0
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 7 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $         XPOS(NWEIGHT), YPOS(NWEIGHT),
     $         FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         PART(NWEIGHT)
            DIVX(NWEIGHT)    = 0.0D0
            DIVY(NWEIGHT)    = 0.0D0
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 9 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $         XPOS(NWEIGHT), YPOS(NWEIGHT),
     $         FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         DIVX(NWEIGHT), DIVY(NWEIGHT),
     $         PART(NWEIGHT)
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 11 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $         XPOS(NWEIGHT), YPOS(NWEIGHT),
     $         FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         DIVX(NWEIGHT), DIVY(NWEIGHT),
     $         CORX(NWEIGHT), CORY(NWEIGHT),
     $         PART(NWEIGHT)

         ELSE
            WRITE(LUNOUT,*) 'SOBP WRONG NO OF COLUMNS',NCOLUMNS
            RETURN
         END IF

      ENDDO


 10   CONTINUE
*     fix index
	  NWEIGHT = NWEIGHT - 1
      WRITE(LUNOUT,*) 'SOBP SOURCE beamlets found:', NWEIGHT

      RETURN
      END

!! ==================================================================================================
!! ==================================================================================================

!> @brief
!! Main user source subroutine
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2010      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA20xy:                                *
*                                                                      *
*     Created on 07 January 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on  17-Oct-10    by    Alfredo Ferrari               *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*  Useful materials on how to prepare user source can be found here:   *
*  http://www.fluka.org/content/course/NEA/lectures/UserRoutines.pdf   *
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*

      INCLUDE 'beamcm.inc'
      INCLUDE 'fheavy.inc'
      INCLUDE 'flkstk.inc'
      INCLUDE 'ioiocm.inc'
      INCLUDE 'ltclcm.inc'
      INCLUDE 'paprop.inc'
      INCLUDE 'sourcm.inc'
      INCLUDE 'sumcou.inc'
      INCLUDE 'caslim.inc'
*
*
*     containers to store data from sobp.dat file
      DOUBLE PRECISION ENERGY(65000), DELTAE(65000)
      DOUBLE PRECISION XPOS(65000), YPOS(65000)
      DOUBLE PRECISION FWHMX(65000), FWHMY(65000)
      DOUBLE PRECISION DIVX(65000), DIVY(65000)
      DOUBLE PRECISION CORX(65000), CORY(65000)
      DOUBLE PRECISION PART(65000)
*
      DOUBLE PRECISION CUMW(65000), TOTW
      SAVE CUMW, TOTW

      INTEGER NWEIGHT, NCOLUMNS
      INTEGER CDF_BINSEARCH
      LOGICAL LPNTSRC
      LOGICAL LENMOMPOS
      DOUBLE PRECISION FWHM2SIGMA
      DOUBLE PRECISION DES
      DOUBLE PRECISION XSPOT, YSPOT
      DOUBLE PRECISION SADX, SADY

      INTEGER IPOS
      INTEGER I, NRAN
      INTEGER IONA
      DOUBLE PRECISION RW, ES, RAN

      CHARACTER*256 FNAME, FPATH

      SAVE LPNTSRC
      SAVE ENERGY, DELTAE, XPOS, YPOS
      SAVE FWHMX, FWHMY, DIVX, DIVY, CORX, CORY, PART
      SAVE NWEIGHT

      LOGICAL LFIRST

      SAVE LFIRST
      DATA LFIRST / .TRUE. /
*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***

* set default file name for sobp.dat if not provided by user in SOURCE card
         FNAME = TRIM(SDUSOU)
         IF ( FNAME.EQ. '' ) THEN
             FNAME = 'sobp.dat'
         END IF

         FWHM2SIGMA = 2.0D0 * SQRT( 2.0D0 * LOG(2.0D0))
         WRITE(LUNOUT,*) 'SPOTLIST SOURCE ZPOS fixed to', ZBEAM
         WRITE(LUNOUT,*) 'SPOTLIST USER PARAM LIST', WHASOU
         WRITE(LUNOUT,*) 'SPOTLIST FNAME', FNAME

         FPATH = '../' // TRIM(FNAME)

*        Fluka run happens in a temporary directory,
*        created in same level as input file
*        sobp.dat is not copied there,
*        We reach one level up to get it via ../sobp.dat
         CALL READSOBP ( TRIM(FPATH), ENERGY, DELTAE,
     $            XPOS, YPOS, FWHMX, FWHMY,
     $            DIVX, DIVY, CORX, CORY, PART, NCOLUMNS, NWEIGHT )

*        In case of problem with reading sobp.dat file
         IF ( (NCOLUMNS .LE. ZERZER) .OR. (NWEIGHT .LE. ZERZER)) THEN
            NOMORE = 1
            RETURN
         ENDIF

*        Build cumulative weights only over strictly positive weights
*        and compact beamlet arrays in-place to exclude zero-weight rows.
         IPOS = 0
         TOTW = 0.0D0
         DO I = 1, NWEIGHT
            IF (PART(I) .LT. 0.0D0) PART(I) = 0.0D0
            IF (PART(I) .GT. 0.0D0) THEN
               IPOS = IPOS + 1
*              Compact all parameter arrays to keep only positive-weight rows
               ENERGY(IPOS) = ENERGY(I)
               DELTAE(IPOS) = DELTAE(I)
               XPOS(IPOS)   = XPOS(I)
               YPOS(IPOS)   = YPOS(I)
               FWHMX(IPOS)  = FWHMX(I)
               FWHMY(IPOS)  = FWHMY(I)
               DIVX(IPOS)    = DIVX(I)
               DIVY(IPOS)    = DIVY(I)
               CORX(IPOS)   = CORX(I)
               CORY(IPOS)   = CORY(I)
               PART(IPOS)   = PART(I)
               TOTW = TOTW + PART(IPOS)
               CUMW(IPOS) = TOTW
            END IF
         END DO
*        Update NWEIGHT to the number of strictly positive-weight rows
         NWEIGHT = IPOS

         IF (TOTW .LE. 0.0D0) THEN
            WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: total weight <= 0'
            NOMORE = 1
            RETURN
         END IF

*        First parameter in SOURCE enables virtual-source geometry.
*        If the number is present and non-zero we expect point-like source
         SADX = WHASOU(3)
         SADY = WHASOU(4)
         IF ( WHASOU(1) .NE. 0.0D0 ) THEN
            WRITE(LUNOUT,*) 'SOBP POINT-LIKE VIRTUAL SOURCE'
            IF ( SADX .LE. 0.0D0 .OR. SADY .LE. 0.0D0 ) THEN
               WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: point-like source'
               WRITE(LUNOUT,*) 'requires positive SADx and SADy, got',
     &            SADX, SADY
               NOMORE = 6
               RETURN
            END IF
            LPNTSRC = .TRUE.
         ELSE
            WRITE(LUNOUT,*) 'SOBP PARALLEL VIRTUAL SOURCE'
            LPNTSRC = .FALSE.
         END IF

      END IF
*  |  End of first call initializations.
*  +-------------------------------------------------------------------*

*  +-------------------------------------------------------------------*
*  |  sample beamlet collection :


*     Get a 64-bit pseudo random number in the interval [0.D+00,1.D+00),
*     1 being not included
      RAN = FLRNDM(111)
      RW  = RAN * TOTW

      NRAN = CDF_BINSEARCH( CUMW, NWEIGHT, RW )

*     check if random number selected properly
      IF ((NRAN .GT. NWEIGHT) .OR. (NRAN .LT. 1)) THEN
         WRITE(LUNOUT,*) 'SOBP SOURCE ERROR. NRAN, RAN:', NRAN, RAN
         NOMORE = 3
         RETURN
      END IF

*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
* Npflka is the stack counter: of course any time source is called it
* must be =0
      NPFLKA = NPFLKA + 1
* Wt is the weight of the particle
      WTFLK (NPFLKA) = 1.0D0   ! particles were already sampled per weight.
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
* Particle type (1=proton.....). Ijbeam is the type set by the BEAM
* card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         IONA   = IPROA
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         IONA   = IPROA
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         IONA  = IBARCH(IONID)
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
      END IF
*  |
*  +-------------------------------------------------------------------*
* From this point .....
* Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
* User dependent flag:
      LOUSE  (NPFLKA) = 0
*  No channeling:
      LCHFLK (NPFLKA) = .FALSE.
      DCHFLK (NPFLKA) = ZERZER
* User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
* User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
* Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
* ... to this point: don't change anything
* Particle age (s)
      AGESTK (NPFLKA) = +ZERZER
      AKNSHR (NPFLKA) = -TWOTWO
****************************************************************



*  +-------------------------------------------------------------------*
*  |  Particle momentum and energy

*      ENERGY and DELTAE from sobp.dat are kinetic energy and sigma
*      in GeV/nucleon. Convert them to total kinetic energy in GeV
*      by multiplying by the integer mass number A, not by the
*      physical ion mass in amu.
        IF (IONA .LE. 0) THEN
           WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: invalid mass number',
     &        IONA, ' for particle ', IONID
           NOMORE = 5
           RETURN
        END IF
        ES   = ENERGY(NRAN)  * DBLE(IONA)
        DES  = DELTAE(NRAN)  * DBLE(IONA)
        IF (DES .LT. 0.0D0) DES = 0.0D0

*  ....................................................................................
*      Basic equation which relates particle energy and momentum is:
*
*         p^2 c^2 + m0^2 c^4 = E^2                 (1)
*
*      Total energy E can be written as the sum of kinetic energy T and rest energy:
*
*         p^2 c^2 + m0^2 c^4 = (T + m0 c^2)^2      (2)
*
*      By doing simple math we can extract from this formula following relations:
*
*         T = sqrt(p^2 c^2 + m0^2 c^4) - m0 c^2    (3)
*         p c = sqrt( E (E + 2 m0 c^2) )           (4)
*
*      We use following variables below:
*        @   PMOFLK - particle momentum, in GeV/c
*        @   AM(IONID) - particle rest energy, in GeV
*        @   ES - total particle kinetic energy, in GeV
*  ....................................................................................
*
*      There is always a small chance that randomly sampled energy or momentum
*      will get negative (i.e. mean energy 10 MeV with energy spread 8 MeV).
*      We will sample gaussian distribution in a loop until positive value is obtained
*
       LENMOMPOS = .FALSE.
       DO WHILE ( LENMOMPOS .EQV. .FALSE. )

          IF ( ES .LE. ZERZER ) THEN
             WRITE(LUNOUT,*) 'SOBP NEGATIVE EN:', ES
             NOMORE = 4
             RETURN
          END IF

*         Lets get a normally distributed random number RGAUSS
          CALL FLNRRN(RGAUSS)

          IF ( NCOLUMNS .GE. 7 ) THEN

*            In case sobp.dat file has 7 columns, we expect than
*            energy spread will be there, set as standard deviation
*            First lets sample gaussian energy distribution.
*            Mean energy is set to ES, value from sobp.dat file
*            Standard deviation is also taken from sobp.dat file
             TKEFLK (NPFLKA) = ES + DES * RGAUSS

*            Exit the loop if kinetic energy is positive (momentum will also be positive)
             IF ( TKEFLK (NPFLKA) .GT. ZERZER ) THEN
                LENMOMPOS = .TRUE.
             ELSE
                WRITE(LUNOUT,*) 'KINETIC ENERGY:', TKEFLK (NPFLKA)
                WRITE(LUNOUT,*) 'RESAMPLING TO GET POSITIVE NUMBER'
             END IF

*            Momentum of the particle, according to eq (4)
             PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA)*
     &        ( TKEFLK (NPFLKA) + TWOTWO * AM (IONID) ))

          ELSE

*            In case energy spread is not provided in sobp.dat file,
*            we can use a momentum spread from BEAM input card (which by default is set to 0).
*            Mean momentum is calculated from kinetic energy of the spot, using eq (4)
*            Standard deviation is obtained from FWHM value set by user in the input card (DPBEAM)
             PMOFLK (NPFLKA) = SQRT ( ES *
     &          ( ES + TWOTWO * AM (IONID) ))
     &          + DPBEAM*RGAUSS/FWHM2SIGMA

*
*            Exit the loop if momentum is positive (kinetic energy will also be positive)
             IF ( PMOFLK (NPFLKA) .GT. ZERZER ) THEN
                LENMOMPOS = .TRUE.
             ELSE
                WRITE(LUNOUT,*) 'MOMENTUM:', PMOFLK (NPFLKA)
                WRITE(LUNOUT,*) 'RESAMPLING TO GET POSITIVE NUMBER'
             END IF

*            Kinetic energy of the particle (GeV), according to eq (3)
             TKEFLK (NPFLKA) = SQRT(PMOFLK(NPFLKA)**2 + AM(IONID)**2)
     &         -AM(IONID)
*

          END IF

       END DO

*      debugging printouts only if requested by user
       IF ( WHASOU(2) .NE. 0.0D0 ) THEN
          WRITE(LUNOUT,*) 'SOBP E:', 1.D3*ES / DBLE(IONA),
     &       'MeV/nucleon'
          WRITE(LUNOUT,*) 'SOBP Z:', ZBEAM, 'cm'
          WRITE(LUNOUT,*) 'SOBP SOURCE EKIN', TKEFLK(NPFLKA)
          WRITE(LUNOUT,*) 'SOBP SOURCE MOMENTUM', PMOFLK(NPFLKA)
       ENDIF

*  +-------------------------------------------------------------------*
*  |  Polarization cosines (TXPOL=-2 flag for "no polarization"):

      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER

C     ------------------------------------------------------------------
C     Reference-plane spot interpretation and SAD source geometry.
C
C     The SOBP input file gives XSPOT and YSPOT as SHIELD-HIT spot
C     coordinates at the reference plane, taken here to be z = 0 cm.
C     They are not used directly as FLUKA particle birth coordinates.
C
C     In point-source/SAD mode, the particle birth position is instead
C     back-projected from the reference-plane spot coordinate to the
C     FLUKA source plane at ZBEAM. The corresponding angular tilt is then
C     applied by APPLY_SAD_TILT below.
C
C     This convention was the source-model correction that made the
C     REFPLANE_MAIN absolute dose and fluence agree with SHIELD-HIT.
C     ------------------------------------------------------------------

*  *  +-------------------------------------------------------------------*
*  |  Particle coordinates + phase space (local beam frame)
      XSPOT = XPOS(NRAN)
      YSPOT = YPOS(NRAN)

*     Interpret XSPOT/YSPOT as SHIELD-HIT reference-plane spot
*     coordinates [cm], not as particle birth coordinates at BEAMPOS.
*     For point-like virtual source/SAD mode, back-project the birth
*     point to the FLUKA source plane ZBEAM [cm]. This makes the
*     SAD-steered central ray pass through the requested spot coordinate
*     near the SHIELD-HIT reference plane z = 0 cm.
      XBEAM = XSPOT
      YBEAM = YSPOT
      IF ( LPNTSRC .AND. SADX .GT. 0.0D0 ) THEN
         XBEAM = XSPOT - (XSPOT / SADX) * (0.0D0 - ZBEAM)
      ENDIF
      IF ( LPNTSRC .AND. SADY .GT. 0.0D0 ) THEN
         YBEAM = YSPOT - (YSPOT / SADY) * (0.0D0 - ZBEAM)
      ENDIF

*     Sample (X,Y,TX,TY) in the local frame (mean angles = 0 here)
      CALL SAMPLE_PHASESPACE(
     &   XBEAM, YBEAM, FWHMX(NRAN), FWHMY(NRAN),
     &   DIVX(NRAN),  DIVY(NRAN),  CORX(NRAN), CORY(NRAN),
     &   FWHM2SIGMA,
     &   XFLK(NPFLKA), YFLK(NPFLKA),
     &   TXFLK(NPFLKA), TYFLK(NPFLKA) )


*     Optional scanning steering (flag-controlled)
      IF ( LPNTSRC ) THEN
         CALL APPLY_SAD_TILT(
     &      XSPOT, YSPOT,
     &      SADX, SADY,
     &      TXFLK(NPFLKA), TYFLK(NPFLKA) )
      ENDIF

*     Renormalize direction
      IF (TXFLK(NPFLKA)**2 + TYFLK(NPFLKA)**2 .GE. 1.0D0) THEN
         TXFLK(NPFLKA) = 0.0D0
         TYFLK(NPFLKA) = 0.0D0
      ENDIF
      TZFLK(NPFLKA) = SQRT(1.0D0
     &   - TXFLK(NPFLKA)**2 - TYFLK(NPFLKA)**2 )

      ZFLK(NPFLKA) = ZBEAM


*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER

*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN
*=== End of subroutine Source =========================================*
      END


      SUBROUTINE SAMPLE_PHASESPACE(
     &   X0, Y0, FWHMX, FWHMY, DIVX, DIVY, CORX, CORY,
     &   FWHM2SIGMA,
     &   X, Y, TX, TY )

      IMPLICIT NONE
      DOUBLE PRECISION X0, Y0, FWHMX, FWHMY, DIVX, DIVY, CORX, CORY
      DOUBLE PRECISION FWHM2SIGMA
      DOUBLE PRECISION X, Y, TX, TY

      DOUBLE PRECISION SIGX, SIGY, SIGTX, SIGTY
      DOUBLE PRECISION RHOX, RHOY, COVX, COVY
      DOUBLE PRECISION G1, G2, G3, G4
      DOUBLE PRECISION DX, DY, DTX, DTY
      DOUBLE PRECISION VARX, VARY

*     Convert spot size FWHM -> sigma (cm)
      SIGX = 0.0D0
      SIGY = 0.0D0
      IF (FWHM2SIGMA .GT. 0.0D0) THEN
         SIGX = FWHMX / FWHM2SIGMA
         SIGY = FWHMY / FWHM2SIGMA
      ENDIF
      IF (SIGX .LT. 0.0D0) SIGX = -SIGX
      IF (SIGY .LT. 0.0D0) SIGY = -SIGY

*     Angular spreads: DIVX/DIVY are in mrad -> convert to rad
      SIGTX = DIVX * 1.0D-3
      SIGTY = DIVY * 1.0D-3
      IF (SIGTX .LT. 0.0D0) SIGTX = -SIGTX
      IF (SIGTY .LT. 0.0D0) SIGTY = -SIGTY

*     CORX/CORY are correlation coefficients rho in [-1,1]
      RHOX = CORX
      RHOY = CORY

*     Clamp rho for numerical safety
      IF (RHOX .GT.  0.999999D0) RHOX =  0.999999D0
      IF (RHOX .LT. -0.999999D0) RHOX = -0.999999D0
      IF (RHOY .GT.  0.999999D0) RHOY =  0.999999D0
      IF (RHOY .LT. -0.999999D0) RHOY = -0.999999D0

*     Convert correlation -> covariance (cm*rad)
      COVX = RHOX * SIGX * SIGTX
      COVY = RHOY * SIGY * SIGTY

*     Draw 4 independent standard normals
      CALL FLNRR2(G1, G2)
      CALL FLNRR2(G3, G4)

*     X-plane: correlated (DX, DTX)
      DX = SIGX * G1

      IF (SIGX .GT. 0.0D0) THEN
         VARX = SIGTX*SIGTX - (COVX*COVX)/(SIGX*SIGX)
         IF (VARX .LT. 0.0D0) VARX = 0.0D0
         DTX = (COVX / SIGX) * G1 + SQRT(VARX) * G2
      ELSE
*        if SIGX=0, correlation is meaningless -> just angle spread
         DTX = SIGTX * G2
      ENDIF

*     Y-plane: correlated (DY, DTY)
      DY = SIGY * G3

      IF (SIGY .GT. 0.0D0) THEN
         VARY = SIGTY*SIGTY - (COVY*COVY)/(SIGY*SIGY)
         IF (VARY .LT. 0.0D0) VARY = 0.0D0
         DTY = (COVY / SIGY) * G3 + SQRT(VARY) * G4
      ELSE
         DTY = SIGTY * G4
      ENDIF

*     Apply to particle (local frame, mean angles = 0)
      X  = X0 + DX
      Y  = Y0 + DY
      TX = DTX
      TY = DTY

      RETURN
      END


      INTEGER FUNCTION CDF_BINSEARCH( CUMW, N, X )
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CUMW(N), X
      INTEGER LO, HI, MID

*     For N <= 1, always return the first (and only) bin
      IF (N .LE. 1) THEN
         CDF_BINSEARCH = 1
         RETURN
      ENDIF

*     Find the first index with CUMW(i) > X (strictly greater).
*     This is a binary search over the monotonic cumulative weights.
      LO = 1
      HI = N
      DO WHILE (LO .LT. HI)
         MID = (LO + HI) / 2
         IF (CUMW(MID) .GT. X) THEN
            HI = MID
         ELSE
            LO = MID + 1
         ENDIF
      END DO

      CDF_BINSEARCH = LO
      RETURN
      END


      SUBROUTINE APPLY_SAD_TILT(XC, YC, SADX, SADY, TX, TY)
      IMPLICIT NONE
      DOUBLE PRECISION XC, YC, SADX, SADY, TX, TY

*     Apply SHIELD-HIT-compatible SAD angular convention.
*     XC and YC are the sampled reference-plane spot coordinates [cm].
*     SADX and SADY are the virtual source-axis distances [cm].
*     Positive XC gives positive TX, i.e. the phase-space envelope spreads
*     from the virtual source through the requested reference-plane spot.
*     This sign, combined with the XBEAM/YBEAM back-projection above,
*     reproduced the REFPLANE_MAIN 10M absolute dose/fluence benchmark.
      IF (SADX .GT. 0.0D0) TX = TX + XC / SADX
      IF (SADY .GT. 0.0D0) TY = TY + YC / SADY

      RETURN
      END

*======================================================================
* LET-moment FLUSCW routine
*======================================================================

*                                                                      *
*=== fluscw ===========================================================*
*                                                                      *
      DOUBLE PRECISION FUNCTION FLUSCW ( IJ    , PLA   , TXX   , TYY   ,
     &                                   TZZ   , WEE   , XX    , YY    ,
     &                                   ZZ    , NREG  , IOLREG, LLO   ,
     &                                   NSURF )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
      INCLUDE 'scohlp.inc'
      INCLUDE 'usrbin.inc'
      INCLUDE 'flkmat.inc'
      INCLUDE 'trackr.inc'
      INCLUDE 'paprop.inc'
      INCLUDE 'fheavy.inc'
      

      DOUBLE PRECISION GETLET
      DOUBLE PRECISION EKIN, LETW, SUMT, SUMD
      INTEGER MATLET, IHEAV, II
      CHARACTER*8 SCONAM


      FLUSCW = ONEONE
      LSCZER = .FALSE.
      SCONAM = TRIM(ADJUSTL(TITUSB(JSCRNG)))

C     ------------------------------------------------------------------
C     Light-fragment LET weighting branches for FLUSCW.
C
C     These branches score LET-weighted fluence contributions for light
C     charged fragments transported by FLUKA:
C
C        IJ = -3    deuteron, 2H
C        IJ = -4    triton, 3H
C        IJ = -5    helium-3, 3He
C        IJ = -6    helium-4 / alpha, 4He
C
C     Scorer-key convention:
C
C        DFL1 / DFL2    deuteron LET / LET^2
C        TFL1 / TFL2    triton LET / LET^2
C        H3L1 / H3L2    helium-3 LET / LET^2
C        H4L1 / H4L2    helium-4 LET / LET^2
C
C     The L1 scorers return one power of LET [keV/um].
C     The L2 scorers return LET^2 [(keV/um)^2].
C
C     Unlike the Li-6/Li-7 branches below, these light fragments use
C     GETLET directly. Lithium is handled separately because GETLET
C     returns zero for the transported lithium ions in this implementation.
C
C     MATLET is taken from MEDFLK(NREG,1), i.e. the material assigned to
C     the current FLUKA region. The material filter below restricts LET
C     scoring to the thesis phantom/slab materials currently expected in
C     this geometry.
C
C     These branches classify the particle currently being transported.
C     They do not record where the fragment was produced or which parent
C     particle produced it. That would require STUPRF or MDSTCK ancestry
C     tagging.
C     ------------------------------------------------------------------

C     Deuteron LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: DFL1.

      IF ( SCONAM .EQ. 'DFL1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -3 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW
         RETURN
      END IF
C     Deuteron LET^2 weighting branch for DLET numerator.
C     Scorer name first four characters: DFL2.

      IF ( SCONAM .EQ. 'DFL2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -3 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW * LETW
         RETURN
      END IF

C     Triton LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: TFL1.

      IF ( SCONAM .EQ. 'TFL1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -4 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW
         RETURN
      END IF
C     Triton LET^2 weighting branch for DLET numerator.
C     Scorer name first four characters: TFL2.

      IF ( SCONAM .EQ. 'TFL2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -4 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW * LETW
         RETURN
      END IF

C     Helium-3 LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: H3L1.

      IF ( SCONAM .EQ. 'H3L1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -5 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW
         RETURN
      END IF
C     Helium-3 LET^2 weighting branch for H3LET numerator.
C     Scorer name first four characters: H3L2.

      IF ( SCONAM .EQ. 'H3L2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -5 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW * LETW
         RETURN
      END IF

C     Helium-4 / alpha LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: H4L1.

      IF ( SCONAM .EQ. 'H4L1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -6 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF
         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW
         RETURN
      END IF
C     Helium-4 / alpha LET^2 weighting branch for H4LET numerator.
C     Scorer name first four characters: H4L2.

      IF ( SCONAM .EQ. 'H4L2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -6 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND.
     &        MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         FLUSCW = LETW * LETW
         RETURN
      END IF

C     ------------------------------------------------------------------
C     Li-6 LET weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L6L1'
C
C     Physics meaning:
C        Keep only transported lithium-6 fragments, identified as
C        charge Z = 3 and mass A = 6, and return one power of LET.
C        This contributes to a track-length-weighted LET numerator.
C
C     FLUKA bookkeeping:
C        JTRACK .LT. -6       means current particle is transported as a
C                             heavy ion / nuclear fragment.
C        NPHEAV .GT. 0        means an entry in FHEAVY is available.
C        KHEAVY(NPHEAV)       maps the current heavy fragment to an
C                             isotope table index IHEAV.
C        ICHEAV(IHEAV)        is fragment charge number Z.
C        IBHEAV(IHEAV)        is fragment mass number A.
C
C     LET reconstruction:
C        GETLET is not used for Li here, because the GETLET call returns
C        zero for these transported lithium ions in this implementation.
C        Instead, local LET is reconstructed from TRACKR step data:
C
C           LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm]
C
C        because 1 GeV/cm = 100 keV/um.
C     ------------------------------------------------------------------


      IF ( SCONAM .EQ. 'L6L1' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-6 LET^2 weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L6L2'
C
C     Physics meaning:
C        Keep only transported lithium-6 fragments, identified as
C        charge Z = 3 and mass A = 6, and return LET^2. Together with
C        the L6L1 scorer, this allows reconstruction of a dose-like
C        or LET-weighted mean for Li-6, depending on the post-processing
C        denominator used.
C
C     Unit reconstruction is identical to L6L1:
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm].
C     ------------------------------------------------------------------

      IF ( SCONAM .EQ. 'L6L2' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW * LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF


C     ------------------------------------------------------------------
C     Li-7 LET weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L7L1'
C
C     Physics meaning:
C        Keep only transported lithium-7 fragments, identified as
C        charge Z = 3 and mass A = 7, and return one power of LET.
C        This is the Li-7 analogue of the L6L1 branch.
C
C     FLUKA isotope identification:
C        JTRACK .LT. -6       selects transported heavy ions/fragments.
C        NPHEAV .GT. 0        requires a valid FHEAVY fragment entry.
C        KHEAVY(NPHEAV)       gives the fragment table index IHEAV.
C        ICHEAV(IHEAV) = 3    requires lithium charge Z = 3.
C        IBHEAV(IHEAV) = 7    requires lithium-7 mass A = 7.
C
C     LET reconstruction:
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm].
C     ------------------------------------------------------------------


      IF ( SCONAM .EQ. 'L7L1' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF

C     ------------------------------------------------------------------
C     Li-7 LET^2 weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L7L2'
C
C     Physics meaning:
C        Keep only transported lithium-7 fragments, identified as
C        charge Z = 3 and mass A = 7, and return LET^2. Together with
C        the L7L1 scorer, this allows reconstruction of a Li-7
C        LET-weighted quantity in post-processing.
C
C     Unit reconstruction is identical to L7L1:
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm].
C     ------------------------------------------------------------------

      IF ( SCONAM .EQ. 'L7L2' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW * LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-6 fluence/filter branch for FLUSCW.
C
C     Scorer key:
C        SCONAM = 'LI6_'
C
C     Physics meaning:
C        Keep only transported lithium-6 fragments, identified as
C        charge Z = 3 and mass A = 6. For matching Li-6 tracks this
C        branch returns ONEONE, so the underlying estimator is scored
C        without additional LET weighting.
C
C     Interpretation:
C        This is an isotope-selection filter. It answers "is the current
C        transported heavy fragment Li-6?" It does not determine where
C        the Li-6 fragment was produced or which parent particle produced
C        it. That ancestry information would require production-time
C        tagging with STUPRF or MDSTCK.
C     ------------------------------------------------------------------

      IF ( SCONAM .EQ. 'LI6_' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN
                  FLUSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-7 fluence/filter branch for FLUSCW.
C
C     Scorer key:
C        SCONAM = 'LI7_'
C
C     Physics meaning:
C        Keep only transported lithium-7 fragments, identified as
C        charge Z = 3 and mass A = 7. For matching Li-7 tracks this
C        branch returns ONEONE, so the underlying estimator is scored
C        without additional LET weighting.
C
C     Interpretation:
C        This is an isotope-selection filter during particle transport.
C        It separates Li-7 from other transported heavy fragments, but
C        it does not record the production vertex, parent particle, or
C        nuclear reaction channel. That would require production-time
C        ancestry tagging with STUPRF or MDSTCK.
C     ------------------------------------------------------------------

      IF ( SCONAM .EQ. 'LI7_' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN
                  FLUSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Proton FLUSCW branch for LET and primary-proton fluence scoring.
C
C     Trigger condition:
C        IJ .EQ. 1          current scored particle is a proton.
C        ISCRNG .EQ. 2      FLUSCW is being called for fluence-like
C                           estimators, including track-length USRBIN.
C
C     Generation convention:
C        LTRACK .EQ. 1      source-generation proton, i.e. one of the
C                           original protons sampled by SOURCE.
C        LTRACK .GT. 1      non-primary proton, i.e. a secondary or later
C                           proton created by a discrete interaction.
C
C     Therefore:
C        PHL1, PHL2, PWL1, PWL2
C           score all transported protons, including source protons and
C           secondary/later-generation protons.
C
C        PRI_, PRL1, PRL2, PWR1, PWR2
C           score only source-generation protons because they require
C           LTRACK .EQ. 1.
C
C     Limitation:
C        LTRACK gives generation number, but not the production vertex,
C        parent particle, target nucleus, or reaction channel. Those
C        ancestry details would require production-time tagging with
C        STUPRF or MDSTCK.
C     ------------------------------------------------------------------

      IF ( IJ .EQ. 1 .AND. ISCRNG .EQ. 2 ) THEN
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF
C        Proton scorer-key map inside this branch:
C
C        PHL1:
C           All-proton LET in the local transport material.
C           Returns LET [keV/um].
C
C        PHL2:
C           All-proton LET^2 in the local transport material.
C           Returns LET^2 [(keV/um)^2].
C
C        PWL1:
C           All-proton LET evaluated in water, independent of local
C           material. Uses MATLET = 30.
C
C        PWL2:
C           All-proton LET^2 evaluated in water. Uses MATLET = 30.
C
C        PRI_:
C           Primary/source-generation proton fluence filter.
C           Returns ONEONE only when LTRACK .EQ. 1.
C
C        PRL1, PRL2:
C           Primary/source-generation proton LET and LET^2 in the local
C           transport material. Require LTRACK .EQ. 1.
C
C        PWR1, PWR2:
C           Primary/source-generation proton LET and LET^2 evaluated in
C           water. Require LTRACK .EQ. 1 and use MATLET = 30.

         IF (SCONAM .EQ. 'PHL1') THEN
            MATLET = MEDFLK(NREG,1)
            IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND. MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
               FLUSCW = ZERZER
               RETURN
            END IF
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            FLUSCW = LETW

         ELSE IF (SCONAM .EQ. 'PHL2') THEN
            MATLET = MEDFLK(NREG,1)
            IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND. MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
               FLUSCW = ZERZER
               RETURN
            END IF
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            FLUSCW = LETW * LETW

         ELSE IF (SCONAM .EQ. 'PWL1') THEN
            MATLET = 30
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            FLUSCW = LETW

         ELSE IF (SCONAM .EQ. 'PWL2') THEN
            MATLET = 30
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            FLUSCW = LETW * LETW

         ELSE IF (SCONAM .EQ. 'PRI_') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               FLUSCW = ONEONE
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM .EQ. 'PRL1') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               MATLET = MEDFLK(NREG,1)
               IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND. MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
                  FLUSCW = ZERZER
                  RETURN
               END IF
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               FLUSCW = LETW
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM .EQ. 'PRL2') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               MATLET = MEDFLK(NREG,1)
               IF ( MATLET .NE. 27 .AND. MATLET .NE. 28 .AND. MATLET .NE. 29 .AND. MATLET .NE. 30 ) THEN
                  FLUSCW = ZERZER
                  RETURN
               END IF
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               FLUSCW = LETW * LETW
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM .EQ. 'PWR1') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               MATLET = 30
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               FLUSCW = LETW
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM .EQ. 'PWR2') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               MATLET = 30
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               FLUSCW = LETW * LETW
            ELSE
               FLUSCW = ZERZER
            END IF

         END IF

      END IF

      RETURN
*=== End of function Fluscw ===========================================*
      END
C=======================================================================
C COMSCW dose-scoring isotope filter.
C=======================================================================
C
C     COMSCW is the user weighting function used by dose-like estimators.
C     This is separate from FLUSCW, which handles fluence-like and
C     track-length estimators.
C
C     Therefore Li-6 and Li-7 isotope filtering for DOSE-like USRBIN
C     scorers must be done here. If this filtering were only present in
C     FLUSCW, the Li-specific dose scorers would not be restricted to the
C     intended lithium isotope.
C
C     Current scope:
C
C        LI6_ / LI6_DZ-like scorer names:
C           keep only transported Li-6 fragments with Z = 3 and A = 6.
C
C        LI7_ / LI7_DZ-like scorer names:
C           keep only transported Li-7 fragments with Z = 3 and A = 7.
C
C     This is transport-time isotope filtering. It does not identify the
C     production vertex, parent particle, target nucleus, or reaction
C     channel. Those ancestry details would require STUPRF or MDSTCK.
C=======================================================================
      DOUBLE PRECISION FUNCTION COMSCW ( IJ    , XA    , YA    , ZA    ,
     &                                   MREG  , RULL  , LLO   , ICALL )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
      INCLUDE 'scohlp.inc'
      INCLUDE 'usrbin.inc'
      INCLUDE 'trackr.inc'
      INCLUDE 'fheavy.inc'

      INTEGER IHEAV
      CHARACTER*8 SCONAM

      LSCZER = .FALSE.
      COMSCW = ONEONE
      SCONAM = TRIM(ADJUSTL(TITUSB(JSCRNG)))

C     ------------------------------------------------------------------
C     Li-6 DOSE filter for COMSCW.
C
C     Scorer key:
C        SCONAM = 'LI6_'
C
C     Physics meaning:
C        For Li-6 dose-like USRBIN scorers, reject every transported
C        particle except lithium-6 fragments. A matching Li-6 fragment is
C        identified by charge Z = 3 and mass A = 6.
C
C     COMSCW return value:
C        COMSCW = ONEONE   keep this energy-deposition contribution.
C        COMSCW = ZERZER   reject this contribution for this scorer.
C
C     Note:
C        This filters dose contributions during transport. It does not
C        identify the production site or parent particle of the Li-6.
C     ------------------------------------------------------------------

      IF ( ISCRNG .EQ. 1 .AND. SCONAM .EQ. 'LI6_' ) THEN
         COMSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN
                  COMSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF

C     ------------------------------------------------------------------
C     Li-7 DOSE filter for COMSCW.
C
C     Scorer key:
C        SCONAM = 'LI7_'
C
C     Physics meaning:
C        For Li-7 dose-like USRBIN scorers, reject every transported
C        particle except lithium-7 fragments. A matching Li-7 fragment is
C        identified by charge Z = 3 and mass A = 7.
C
C     COMSCW return value:
C        COMSCW = ONEONE   keep this energy-deposition contribution.
C        COMSCW = ZERZER   reject this contribution for this scorer.
C
C     Note:
C        This filters dose contributions during transport. It does not
C        identify the production site or parent particle of the Li-7.
C     ------------------------------------------------------------------

      IF ( ISCRNG .EQ. 1 .AND. SCONAM .EQ. 'LI7_' ) THEN
         COMSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN
                  COMSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF

      RETURN
*=== End of function Comscw ===========================================*
      END