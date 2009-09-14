export FFLAGS="-arch x86_64"
export LDFLAGS="-arch x86_64 -undefined dynamic_lookup -dynamiclib"
export MACOSX_DEPLOYMENT_TARGET=10.6

maker_dir := StRoot/StSpinPool/StMCAsymMaker

objects := _aac06.so _bb.so _ctq5par.so _dns.so _dssv.so _grsv.so \
	_grsv2000.so _grv.so _lss2006.so _polnlo.so

all: $(objects)

clean:
	rm *.so

_aac06.so: aac06.pyf $(maker_dir)/aac06.F
	f2py -c $?;

_bb.so: bb.pyf $(maker_dir)/ppdf_bb.F
	f2py -c $?;

_ctq5par.so: ctq5par.pyf $(maker_dir)/Ctq5Par.F
	f2py -c $?;

_dns.so: dns.pyf $(maker_dir)/readgrid_dns.F
	f2py -c $?;

_dssv.so: dssv.pyf $(maker_dir)/DSSV.F
	f2py -c $?;

_grsv.so: grsv.pyf $(maker_dir)/grsv.F
	f2py -c $?;

_grsv2000.so: grsv2000.pyf $(maker_dir)/grsv2000std_val_LO_NLO.F
	f2py -c $?;

_grv.so: grv.pyf $(maker_dir)/grv.F
	f2py -c $?;

_lss2006.so: lss2006.pyf $(maker_dir)/LSS2006pdf_g1.F
	f2py -c $?;

_polnlo.so: polnlo.pyf $(maker_dir)/polnlo.F
	f2py -c $?;
