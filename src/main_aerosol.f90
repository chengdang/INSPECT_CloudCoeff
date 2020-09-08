!-------------------------------------------------------
!
! main_aerosol.f90
!
! Description:
! ============
!	Simple test program to inspect the CRTM AerosolCoeff
!	binary coefficient files that contain the scattering
! properties for aerosols.
! For a definition of the AerosolCoeff_type TYPE, please
! look at AerosolCoeff_Define.f90.
!
! SIDE-EFFECTS:
! =============
! Overwrites all existing ASCII files in the folder output_aerosol.
!
!	Record of Revisions:
! ====================
!
!	Date: 	    Author:        Description:
! =====       =======        ============
! 2020-07-30  C. Dang        Aerosol coefficient inspect code
!                            Code description and format follow main.f90
!                            by Patrick Stegmann
!
! Copyright © 2018 Patrick Stegmann
!
! This file is part of INSPECT_AerosolCoeff.

! This program is distributed in the hope that it will
! be useful,
! but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
! PURPOSE.  See the Apache License for more details.
!
! You should have received a copy of the Apache 2.0
! License along with this program. If not,
! see <https://www.apache.org/licenses/LICENSE-2.0>.
!
!-------------------------------------------------------


!-------------------------------------------------------

PROGRAM inspect_AerosolCoeff
    ! ============================================================================
    ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
    !
    ! Module usage
    USE CRTM_Module
    USE AerosolCoeff_Binary_IO, ONLY: AerosolCoeff_Binary_InquireFile, &
                                      AerosolCoeff_Binary_ReadFile, &
                                      AerosolCoeff_Binary_WriteFile
    USE AerosolCoeff_netCDF_IO, ONLY: AerosolCoeff_netCDF_InquireFile, &
                                      AerosolCoeff_netCDF_ReadFile, &
                                      AerosolCoeff_netCDF_WriteFile
    USE AerosolCoeff_Define,    ONLY: AerosolCoeff_type, &
                                      AerosolCoeff_Associated, &
                                      AerosolCoeff_Destroy, &
                                      AerosolCoeff_Inspect

    ! Disable all implicit typing
    IMPLICIT NONE
    ! ============================================================================

    ! String lengths
    INTEGER, PARAMETER :: ML = 256   ! Error message length
    INTEGER, PARAMETER :: SL = 500  ! Maximum length for path+filenames

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'AerosolCoeff test read'

    ! Arguments
    INTEGER       :: Process_ID
    INTEGER       :: Output_Process_ID

    INTEGER :: err_stat
    INTEGER :: Destroy_Status
    LOGICAL :: noisy
    LOGICAL :: Quiet
    LOGICAL :: selection
    INTEGER :: aerosol_selector, &
               wvl_selector, &
               rad_selector

    ! Local variables
    CHARACTER(ML)   :: msg, pid_msg
    CHARACTER(SL)   :: Bin_Default_AerosolCoeff_File, Bin_AerosolCoeff_File
    CHARACTER(SL)   :: NC_Default_AerosolCoeff_File, NC_AerosolCoeff_File
    CHARACTER(SL)   :: Filename
    CHARACTER(SL)   :: File_Path
    CHARACTER(SL)   :: FP_output

    INTEGER ::   n_Wavelengths   , &
                 n_Radii         , &
                 n_Types         , &
                 n_RH            , &
                 n_Legendre_Terms, &
                 n_Phase_Elements, &
                 Release         , &
                 Version

    INTEGER :: ii, jj ! Iterators


    REAL, DIMENSION(31,10) :: w, g, ke, w_nc, g_nc, ke_nc
    REAL, DIMENSION(31,10,38) :: pcoeff, pcoeff_nc
    TYPE(AerosolCoeff_type), TARGET, SAVE :: AeroC, AeroC_NC


! ---------------------------------
!  Binary LUT
! ---------------------------------
    ! Set up
    WRITE(*,*) "Binary IO"
    err_stat = SUCCESS
    Bin_Default_AerosolCoeff_File    = 'AerosolCoeff.bin'
    Filename  = Bin_Default_AerosolCoeff_File
    File_Path  = '../fix/'
    FP_output = './output_aerosol/Binary'

    ! ...Add the file path
    Bin_AerosolCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(Filename)

    ! Inquire the content of the AerosolCoeff data file
    err_stat = AerosolCoeff_Binary_InquireFile(&
                 Bin_AerosolCoeff_File, &
                 n_Wavelengths     = n_Wavelengths   , &
                 n_Radii          = n_Radii          , &
                 n_Types          = n_Types          , &
                 n_RH             = n_RH             , &
                 n_Legendre_Terms = n_Legendre_Terms , &
                 n_Phase_Elements = n_Phase_Elements , &
                 Release          = Release          , &
                 Version          = Version            )

    ! Write General Info about the AerosolCoeff.bin to ASCII file.
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//TRIM('General_Info.txt'))
    WRITE(1,*) "n_Wavelengths = ", n_Wavelengths, NEW_LINE('A'),&
               "n_Radii = ", n_Radii, NEW_LINE('A'),&
               "n_Types = ", n_Types, NEW_LINE('A'),&
               "n_RH = ", n_RH, NEW_LINE('A'),&
               "n_Legendre_Terms = ", n_Legendre_Terms, NEW_LINE('A'),&
               "n_Phase_Elements = ", n_Phase_Elements, NEW_LINE('A'),&
               "Release: ", Release, NEW_LINE('A'),&
               "Version: ", Version
    CLOSE(1)


    ! Read the CloudCoeff data file
    err_stat = AerosolCoeff_Binary_ReadFile( &
               Bin_AerosolCoeff_File, &
               AeroC, &
               Quiet = .NOT. noisy )
    IF ( err_stat /= SUCCESS ) THEN
      WRITE( msg,'("Error reading AerosolCoeff file ",a)') TRIM(Bin_AerosolCoeff_File)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    ELSE
      WRITE( msg,'("Success reading AerosolCoeff file ",a)') TRIM(Bin_AerosolCoeff_File)
    END IF


    !-------------------------------------------------------
    !
    ! CRTM Bulk Scattering Properties:
    !
    ! Abreviations:
    !   Reff: Effective radius
    !   ke:   Extinction coefficient
    !   w:    Single scatter albedo
    !   g:    Asymmetry parameter
    !
    !-------------------------------------------------------

    !-------------------------------------------------------
    ! Query Particle type as defined by CRTM
    !-------------------------------------------------------

    ! Select aerosol type
    WRITE(*,*) "Please select the aerosol type by number:"
    WRITE(*,*) "1: Dust"
    WRITE(*,*) "2: Sea salt-SSAM"
    WRITE(*,*) "3: Sea salt-SSCM1"
    WRITE(*,*) "4: Sea salt-SSCM2"
    WRITE(*,*) "5: Sea salt-SSCM3"
    WRITE(*,*) "6: Organic carbon"
    WRITE(*,*) "7: Black carbon"
    WRITE(*,*) "8: Sulfate"
    READ(*,*)  aerosol_selector
    WRITE(*,*) "Selected aerosol type is: ", AeroC%Type_Name(aerosol_selector)

    ! Select wavelength
    WRITE(*,*) "To retrieve the corresponding Legendre coefficients,"
    WRITE(*,*) "Please select the wavelength as an integer between 1 and ", n_Wavelengths
    READ(*,*) wvl_selector
    WRITE(*,*) "Selected wavelength is [microns]: ", AeroC%Wavelength(wvl_selector)

    ! Select effective radius
    WRITE(*,*) "and the effective radius as an integer between 1 and ", n_Radii
    READ(*,*) rad_selector
    WRITE(*,*) "Selected radius is [microns]: ", AeroC%Reff(rad_selector, aerosol_selector)


    !-------------------------------------------------------
    !
    ! Output AerosolCoeff data to file:
    !
    !-------------------------------------------------------

    !-------------------------------------------------------
    ! Wavelength
    !-------------------------------------------------------

    ! Write selected index/properties for NetCDF file use
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//"netCDF_information.txt")
      WRITE(1,*) aerosol_selector, wvl_selector, rad_selector
      ! WRITE(1,*) wvl_selector, AeroC%Wavelength(wvl_selector)
      ! WRITE(1,*) rad_selector, AeroC%Reff(rad_selector, aerosol_selector)
    CLOSE(1)

    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Wavelength.txt')
    !WRITE(1,*) 'CRTM aerosol scattering property wavelength in units of [Microns (um)].'
    DO ii = 1,n_Wavelengths
        WRITE(1,*) AeroC%Wavelength(ii)
    END DO
    CLOSE(1)

    !-------------------------------------------------------
    ! Scattering Particle Effective Radii
    !-------------------------------------------------------
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Radii.txt')
    DO ii = 1,n_Radii
        WRITE(1,*) AeroC%Reff(ii, aerosol_selector)
    END DO
    CLOSE(1)


    ! Write Legendre coefficients to file.
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//"Legendre_terms.txt")
    DO ii = 0,n_Legendre_Terms
        WRITE(1,*) AeroC%pcoeff(wvl_selector, rad_selector, aerosol_selector, ii, 1)
    END DO
    CLOSE(1)

    ! Write extinction coefficients
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Extinction_Coefficients.txt')
      WRITE(1,*) AeroC%ke(:,:,aerosol_selector)
    CLOSE(1)

    ! Write single scattering albedo
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'SingleScatAlbedo.txt')
      WRITE(1,*) AeroC%w(:,:,aerosol_selector)
    CLOSE(1)

    ! Write asymmetry facyor
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Asymmery_factor.txt')
      WRITE(1,*) AeroC%g(:,:,aerosol_selector)
    CLOSE(1)



    ! ---------------------------------
    !  netCDF LUT
    ! ---------------------------------
    ! Set up
    WRITE(*,*) "NetCDF IO"
    err_stat = SUCCESS
    NC_Default_AerosolCoeff_File    = 'AerosolCoeff.nc'
    Filename  = NC_Default_AerosolCoeff_File
    File_Path  = '../fix/'
    FP_output = './output_aerosol/NetCDF'

    ! ...Add the file path
    NC_AerosolCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(Filename)

    ! Inquire the content of the AerosolCoeff data file
    err_stat = AerosolCoeff_netCDF_InquireFile(&
                 NC_AerosolCoeff_File, &
                 n_Wavelengths     = n_Wavelengths   , &
                 n_Radii          = n_Radii          , &
                 n_Types          = n_Types          , &
                 n_RH             = n_RH             , &
                 n_Legendre_Terms = n_Legendre_Terms , &
                 n_Phase_Elements = n_Phase_Elements , &
                 Release          = Release          , &
                 Version          = Version            )

    ! Write General Info about the AerosolCoeff.bin to ASCII file.
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//TRIM('General_Info.txt'))
    WRITE(1,*) "n_Wavelengths = ", n_Wavelengths, NEW_LINE('A'),&
               "n_Radii = ", n_Radii, NEW_LINE('A'),&
               "n_Types = ", n_Types, NEW_LINE('A'),&
               "n_RH = ", n_RH, NEW_LINE('A'),&
               "n_Legendre_Terms = ", n_Legendre_Terms, NEW_LINE('A'),&
               "n_Phase_Elements = ", n_Phase_Elements, NEW_LINE('A'),&
               "Release: ", Release, NEW_LINE('A'),&
               "Version: ", Version
    CLOSE(1)


    ! Read the CloudCoeff data file
    err_stat = AerosolCoeff_netCDF_ReadFile( &
               NC_AerosolCoeff_File, &
               AeroC_NC, &
               Quiet = .NOT. noisy )
    IF ( err_stat /= SUCCESS ) THEN
      WRITE( msg,'("Error reading AerosolCoeff file ",a)') TRIM(NC_AerosolCoeff_File)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    ELSE
      WRITE( msg,'("Success reading AerosolCoeff file ",a)') TRIM(NC_AerosolCoeff_File)
    END IF


    !-------------------------------------------------------
    !
    ! CRTM Bulk Scattering Properties:
    !
    ! Abreviations:
    !   Reff: Effective radius
    !   ke:   Extinction coefficient
    !   w:    Single scatter albedo
    !   g:    Asymmetry parameter
    !
    !-------------------------------------------------------

    !-------------------------------------------------------
    ! Query Particle type as defined by CRTM
    !-------------------------------------------------------

    ! Select aerosol type
    WRITE(*,*) "NetCDF IO"
    WRITE(*,*) "Please select the aerosol type by number:"
    WRITE(*,*) "1: Dust"
    WRITE(*,*) "2: Sea salt-SSAM"
    WRITE(*,*) "3: Sea salt-SSCM1"
    WRITE(*,*) "4: Sea salt-SSCM2"
    WRITE(*,*) "5: Sea salt-SSCM3"
    WRITE(*,*) "6: Organic carbon"
    WRITE(*,*) "7: Black carbon"
    WRITE(*,*) "8: Sulfate"
    READ(*,*)  aerosol_selector
    WRITE(*,*) "Selected aerosol type is: ", AeroC_NC%Type_Name(aerosol_selector)

    ! Select wavelength
    WRITE(*,*) "To retrieve the corresponding Legendre coefficients,"
    WRITE(*,*) "Please select the wavelength as an integer between 1 and ", n_Wavelengths
    READ(*,*) wvl_selector
    WRITE(*,*) "Selected wavelength is [microns]: ", AeroC_NC%Wavelength(wvl_selector)

    ! Select effective radius
    WRITE(*,*) "and the effective radius as an integer between 1 and ", n_Radii
    READ(*,*) rad_selector
    WRITE(*,*) "Selected radius is [microns]: ", AeroC_NC%Reff(rad_selector, aerosol_selector)


    !-------------------------------------------------------
    !
    ! Output AerosolCoeff data to file:
    !
    !-------------------------------------------------------

    !-------------------------------------------------------
    ! Wavelength
    !-------------------------------------------------------

    ! Write selected index/properties for NetCDF file use
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//"netCDF_information.txt")
      WRITE(1,*) aerosol_selector, wvl_selector, rad_selector
      ! WRITE(1,*) wvl_selector, AeroC%Wavelength(wvl_selector)
      ! WRITE(1,*) rad_selector, AeroC%Reff(rad_selector, aerosol_selector)
    CLOSE(1)

    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Wavelength.txt')
    !WRITE(1,*) 'CRTM aerosol scattering property wavelength in units of [Microns (um)].'
    DO ii = 1,n_Wavelengths
        WRITE(1,*) AeroC_NC%Wavelength(ii)
    END DO
    CLOSE(1)

    !-------------------------------------------------------
    ! Scattering Particle Effective Radii
    !-------------------------------------------------------
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Radii.txt')
    DO ii = 1,n_Radii
        WRITE(1,*) AeroC_NC%Reff(ii, aerosol_selector)
    END DO
    CLOSE(1)


    ! Write Legendre coefficients to file.
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//"Legendre_terms.txt")
    DO ii = 0,n_Legendre_Terms
        WRITE(1,*) AeroC_NC%pcoeff(wvl_selector, rad_selector, aerosol_selector, ii, 1)
    END DO
    CLOSE(1)

    ! Write extinction coefficients
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Extinction_Coefficients.txt')
      WRITE(1,*) AeroC_NC%ke(:,:,aerosol_selector)
    CLOSE(1)

    ! Write single scattering albedo
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'SingleScatAlbedo.txt')
      WRITE(1,*) AeroC_NC%w(:,:,aerosol_selector)
    CLOSE(1)

    ! Write asymmetry facyor
    OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FP_output))//'Asymmery_factor.txt')
      WRITE(1,*) AeroC_NC%g(:,:,aerosol_selector)
    CLOSE(1)

END PROGRAM inspect_AerosolCoeff
