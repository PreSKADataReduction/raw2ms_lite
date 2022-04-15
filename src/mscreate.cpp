// MSCreate.cc: Implementation for creating a MeasurementSet
//
//  Copyright (C) 2005
//  ASTRON (Netherlands Foundation for Research in Astronomy)
//  P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, seg@astron.nl
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  $Id: MSCreate.cc 12842 2009-03-09 08:11:29Z diepen $


#include "../include/mscreate.hpp"
using namespace ulastai;
using namespace casacore;


raw_data_source::raw_data_source(int n)
  :nbaselines(n)
{}

int raw_data_source::num_of_baselines()const
{
  return nbaselines;
}

bool raw_data_source::fetch_one()
{
  return do_fetch_one();
}

std::pair<int,int> raw_data_source::antenna_pair(int bl)const
{
  return do_antenna_pair(bl);
}

Array<Complex> raw_data_source::data(int field,int band,int bl)const
{
  return do_data(field,band,bl);
}

Array<Float> raw_data_source::sigma(int field,int band,int bl)const
{
  return do_sigma(field,band,bl);
}

Array<Bool> raw_data_source::flags(int field,int band,int bl)const
{
  return do_flags(field,band,bl);
}

double raw_data_source::time()const
{
  return do_time();
}

mscreate::mscreate (const std::string& ms_name,
		    double start_time, int npol,
		    const Table& ant_tab,
		    const casacore::MPosition& array_pos)
  : its_nbands        (0),
    its_nfields       (0),
    its_nantennas         (0),
    its_npol_per_ant        (npol),
    //itsNrTimes       (0),
    //itsTimeStep      (timeStep),
    its_start_time     (start_time),
    its_end_time       (start_time),
    its_npol         (new Block<Int>),
    its_nchan        (new Block<Int>),
    its_poln         (new Block<Int>),
    its_ant_bl         (0),
    its_array_pos(new casacore::MPosition(array_pos)),
    its_frame         (new MeasFrame(*its_array_pos)),
    its_phase_dir      (new Block<MDirection>()),
    its_ms            (0),
    its_ms_col         (0),
    correct_w(false)
{
  
  // Use the middle antenna as the array position.
  // Setup the frame for the UVW calculations.
  //its_array_pos = new MPosition(antMPos[nantennas/2]);
  //its_frame = new MeasFrame(*its_array_pos);
  // Create the MS.
  
  create_ms (ms_name, ant_tab);
  // Fill the baseline vector for each antenna.
  fill_baselines();
}

mscreate::~mscreate()
{
  if (its_ms != 0) {
    update_times();
  }
  delete its_npol;
  delete its_nchan;
  delete its_poln;
  delete its_ant_bl;
  delete its_array_pos;
  delete its_frame;
  delete its_phase_dir;
  delete its_ms_col;
  delete its_ms;
}

int mscreate::num_of_polarizations() const
{
  return its_ms->polarization().nrow();
}


void mscreate::create_ms (const String& ms_name,
			  //const Block<MPosition>& antPos,
			  const Table& ant_tab)
{
  // Create an integer flag column?
  //IPosition dataShape(2,its_npol_per_ant,itsNrFreq);
  // Get the MS main default table description.
  TableDesc td = MS::requiredTableDesc();
  // Add the data column and its unit.
  MS::addColumnToDesc(td, MS::DATA, 2);
  td.rwColumnDesc(MS::columnName(MS::DATA)).rwKeywordSet().
    define("UNIT","Jy");
  //Vector<String> tsmNames(1);
  //tsmNames[0] = MS::columnName(MS::DATA);
  //td.rwColumnDesc(tsmNames[0]).setShape (dataShape);
  //td.defineHypercolumn("TiledData", 3, tsmNames);
  //tsmNames[0] = MS::columnName(MS::FLAG);
  //td.rwColumnDesc(tsmNames[0]).setShape (dataShape);
  //td.defineHypercolumn("TiledFlag", 3, tsmNames);////
  //if (n_flag_bits <= 0) {
  //  td.defineHypercolumn("TiledFlag", 3, tsmNames);////
  //}
  //else {
  //tsmNames[0] = flag_column;
  //td.rwColumnDesc(flag_column).setShape (dataShape);
  //td.defineHypercolumn("TiledIntFlag", 3, tsmNames);////
  //}
  ////td.defineHypercolumn("TiledFlag", 3, tsmNames);
  //tsmNames[0] = MS::columnName(MS::UVW);
  //td.defineHypercolumn("TiledUVW", 2, tsmNames);
  // Setup the new table.
  // Most columns use the IncrStMan; some use others.
  SetupNewTable new_tab(ms_name, td, Table::New);
  IncrementalStMan incr_std_man("ISMData");
  new_tab.bindAll (incr_std_man);
  StandardStMan    stan_std_man;
  new_tab.bindColumn(MS::columnName(MS::ANTENNA1), stan_std_man);
  new_tab.bindColumn(MS::columnName(MS::ANTENNA2), stan_std_man);
  // Use a TiledColumnStMan for the data, flags and UVW.
  // Store all pol and freq in a single tile.
  // In this way the data appear in separate files that can be mmapped.
  // Flags are stored as bits, so take care each tile has multiple of 8 flags.
  // Create the MS and its subtables.
  // Get access to its columns.
  its_ms = new MeasurementSet(new_tab);
  its_ms_col = new MSMainColumns(*its_ms);
  // Create all subtables.
  // Do this after the creation of optional subtables,
  // so the MS will know about those optional sutables.
  its_ms->createDefaultSubtables (Table::New);
  // Fill various subtables.
  fill_antenna (ant_tab);
  fill_feed();
  fill_processor();
  fill_observation();
  fill_state();
  // Find out which datamanagers contain DATA, FLAG and UVW.
  // Create symlinks for them, but only if no tiling in frequency.
}

void mscreate::set_correct_w(bool b)
{
  correct_w=b;
}

int mscreate::add_band (int nchannels,
		       double ref_freq, double chan_width)
{
  AlwaysAssert (nchannels > 0, AipsError);
  Vector<double> chan_widths(nchannels);
  chan_widths = chan_width;
  Vector<double> chan_freqs(nchannels);
  indgen (chan_freqs, ref_freq - (nchannels-1)*chan_width/2., chan_width);
  return add_band (nchannels, ref_freq, chan_freqs, chan_widths);
}

int mscreate::add_band (int nchannels,
		       double ref_freq, const double* chan_freqs,
		       const double* chan_widths)
{
  AlwaysAssert (nchannels > 0, AipsError);
  IPosition shape(1, nchannels);
  Vector<double> freqs (shape, const_cast<double*>(chan_freqs), SHARE);
  Vector<double> widths(shape, const_cast<double*>(chan_widths), SHARE);
  return add_band (nchannels, ref_freq, freqs, widths);
}

int mscreate::add_band (int nchannels,
		       double ref_freq, const Vector<double>& chan_freqs,
		       const Vector<double>& chan_widths)
{
  std::vector<double> f(nchannels);
  for(int i=0;i<nchannels;++i)
    {
      f[i]=chan_freqs[i];
    }
  ch_freq_vectors.push_back(f);
  int npolarizations=its_npol_per_ant;
  AlwaysAssert (npolarizations==1 || npolarizations==2 || npolarizations==4,
		AipsError);
  AlwaysAssert (nchannels > 0, AipsError);
  //AlwaysAssert (npolarizations==its_npol_per_ant && nchannels==itsNrFreq, AipsError);
  // Find out if this nr of polarizations has already been given.
  Int polnr = -1;
  for (Int i=0; i<its_nbands; i++) {
    if (npolarizations == (*its_npol)[i]) {
      polnr = (*its_poln)[i];
      break;
    }
  }
  // If not, add an entry to the POLARIZATION subtable.
  if (polnr < 0) {
    polnr = add_polarization (npolarizations);
  }

  
  // Add a row to the DATADESCRIPTION subtable.
  MSDataDescription msdd = its_ms->dataDescription();
  MSDataDescColumns msdd_col(msdd);
  uInt rownr = msdd.nrow();
  msdd.addRow();
  msdd_col.spectralWindowId().put (rownr, rownr);
  msdd_col.polarizationId().put (rownr, polnr);
  msdd_col.flagRow().put (rownr, False);
  // Add a row to the SPECTRALWINDOW subtable.
  // Find the total bandwidth from the minimum and maximum.
  Vector<double> start_freqs = chan_freqs - chan_widths/2.;
  Vector<double> end_freqs = chan_freqs + chan_widths/2.;
  double total_bw = max(end_freqs) - min(start_freqs);
  MSSpectralWindow ms_spw = its_ms->spectralWindow();
  MSSpWindowColumns ms_spw_col(ms_spw);
  ms_spw.addRow();
  ms_spw_col.numChan().put (rownr, nchannels);
  ms_spw_col.name().put (rownr, "");
  ms_spw_col.refFrequency().put (rownr, ref_freq);
  ms_spw_col.chanFreq().put (rownr, chan_freqs);
  ms_spw_col.chanWidth().put (rownr, chan_widths);
  ms_spw_col.measFreqRef().put (rownr, MFrequency::TOPO);
  ms_spw_col.effectiveBW().put (rownr, chan_widths);
  ms_spw_col.resolution().put (rownr, chan_widths);
  ms_spw_col.totalBandwidth().put (rownr, total_bw);
  ms_spw_col.netSideband().put (rownr, 0);
  ms_spw_col.ifConvChain().put (rownr, 0);
  ms_spw_col.freqGroup().put (rownr, 0);
  ms_spw_col.freqGroupName().put (rownr, "");
  ms_spw_col.flagRow().put (rownr, False);
  // Now add the band to the internal blocks.
  its_nbands++;
  its_npol->resize (its_nbands);
  its_nchan->resize (its_nbands);
  its_poln->resize (its_nbands);
  (*its_npol)[its_nbands-1] = npolarizations;
  (*its_nchan)[its_nbands-1] = nchannels;
  (*its_poln)[its_nbands-1] = polnr;
  ms_spw.flush();
  msdd.flush();
  return its_nbands-1;
}

int mscreate::add_polarization (int npolarizations)
{
  MSPolarization ms_pol = its_ms->polarization();
  MSPolarizationColumns ms_pol_col(ms_pol);
  uInt rownr = ms_pol.nrow();
  Vector<Int> corr_type(npolarizations);
  if (npolarizations == 1)
  {
    corr_type(0) = Stokes::I;
  }
  else if (npolarizations == 2)
  {
    corr_type(0) = Stokes::XX;
    corr_type(1) = Stokes::YY;
  }
  else if (npolarizations == 4)
  {
    corr_type(1) = Stokes::XY;
    corr_type(2) = Stokes::YX;
    corr_type(3) = Stokes::YY;
  }
  else
  {
    assert(0);
  }
  Matrix<Int> corr_product(2, npolarizations);
  for (Int i=0; i<npolarizations; i++) {
    corr_product(0,i) = Stokes::receptor1(Stokes::type(corr_type(i)));
    corr_product(1,i) = Stokes::receptor2(Stokes::type(corr_type(i)));
  }
  // Fill the columns.
  ms_pol.addRow();
  ms_pol_col.numCorr().put (rownr, npolarizations);
  ms_pol_col.corrType().put (rownr, corr_type);
  ms_pol_col.corrProduct().put (rownr, corr_product);
  ms_pol_col.flagRow().put (rownr, False);
  ms_pol.flush();
  return rownr;
}

int mscreate::add_field (double ra, double dec)
{
  MVDirection radec (Quantity(ra,"rad"), Quantity(dec,"rad"));
  MDirection indir(radec, MDirection::J2000);
  //if (its_phase_dir == 0) {
  //  its_phase_dir = new Block<MDirection>();
  //}
  its_phase_dir->resize (its_nfields+1);
  (*its_phase_dir)[its_nfields] = indir;
  Vector<MDirection> outdir(1);
  outdir(0) = indir;
  // Put the direction into the FIELD subtable.
  {
    MSField ms_field = its_ms->field();
    MSFieldColumns ms_field_col(ms_field);
    uInt rownr = ms_field.nrow();
    ms_field.addRow();
    ms_field_col.name().put (rownr, "BEAM_" + String::toString(rownr));
    ms_field_col.code().put (rownr, "");
    ms_field_col.time().put (rownr, its_start_time);
    ms_field_col.numPoly().put (rownr, 0);
    ms_field_col.delayDirMeasCol().put (rownr, outdir);
    ms_field_col.phaseDirMeasCol().put (rownr, outdir);
    ms_field_col.referenceDirMeasCol().put (rownr, outdir);
    ms_field_col.sourceId().put (rownr, -1);
    ms_field_col.flagRow().put (rownr, False);
  }
  // Put the direction for each antenna into the POINTING subtable.
  {
    MSPointing ms_pointing = its_ms->pointing();
    MSPointingColumns ms_pointing_col(ms_pointing);
    uInt rownr = ms_pointing.nrow();
    ms_pointing.addRow(its_nantennas);
    for (Int i=0; i<its_nantennas; i++) {
      ms_pointing_col.antennaId().put (rownr, i);
      ms_pointing_col.time().put (rownr, its_start_time);
      ms_pointing_col.interval().put (rownr, 0.);
      ms_pointing_col.name().put (rownr, "");
      ms_pointing_col.numPoly().put (rownr, 0);
      ms_pointing_col.timeOrigin().put (rownr, its_start_time);
      ms_pointing_col.directionMeasCol().put (rownr, outdir);
      ms_pointing_col.targetMeasCol().put (rownr, outdir);
      ms_pointing_col.tracking().put (rownr, False);
      rownr++;
    }
  }
  its_nfields++;
  return its_nfields-1;
}

void mscreate::fill_antenna (const Table& ant_tab)
{
  ROArrayColumn<double> pos_col(ant_tab, "POSITION");
  ROArrayColumn<double> offset_col(ant_tab,"OFFSET");
  Array<double> its_ant_pos(pos_col.getColumn());
  Matrix<double> its_offset(Array<double>(offset_col.getColumn()));
  its_nantennas=its_ant_pos.shape()[1];
  Matrix<double> ant_pos(its_ant_pos);
  Block<MPosition> ant_mpos(its_nantennas);
  for(int i=0;i<its_nantennas;++i)
    {
      ant_mpos[i]=MPosition(MVPosition(ant_pos(0,i),ant_pos(1,i),ant_pos(2,i)),
			   MPosition::ITRF);
    }
  //its_array_pos=new MPosition(ant_mpos[its_nantennas/2]);
  //its_frame = new MeasFrame(*its_array_pos);

  ROScalarColumn<String> name_col(ant_tab,"NAME");
  Vector<String> its_ant_name(name_col.getColumn());

  ROScalarColumn<String> station_col(ant_tab,"STATION");
  Vector<String> its_station_name(station_col.getColumn());

  ROScalarColumn<String> mount_col(ant_tab,"MOUNT");
  Vector<String> its_mount(mount_col.getColumn());

  ROScalarColumn<double> diameter_col(ant_tab,"DISH_DIAMETER");
  Vector<double> its_diameter(diameter_col.getColumn());
  
  // Determine constants for the ANTENNA subtable.
  Vector<Double> ant_offset(3);
  ant_offset = 0;
  // Fill the ANTENNA subtable.
  MSAntenna ms_ant = its_ms->antenna();
  MSAntennaColumns ms_ant_col(ms_ant);
  ms_ant.addRow (its_nantennas);
  for (Int i=0; i<its_nantennas; i++) {
    //ms_ant_col.name().put (i, "ST_" + String::toString(i));
    ms_ant_col.name().put (i, its_ant_name[i]);
    ms_ant_col.station().put (i, its_station_name[i]);
    ms_ant_col.type().put (i, "GROUND-BASED");
    ms_ant_col.mount().put (i, its_mount[i]);
    ms_ant_col.positionMeas().put (i, ant_mpos[i]);
    ant_offset[0]=its_offset(0,i);
    ant_offset[1]=its_offset(1,i);
    ant_offset[2]=its_offset(2,i);
    ms_ant_col.offset().put (i, ant_offset);
    ms_ant_col.dishDiameter().put (i, its_diameter[i]);
    ms_ant_col.flagRow().put (i, False);
  }
  ms_ant.flush();
}

void mscreate::fill_feed()
{
  // Determine constants for the FEED subtable.
  Int nRec = 2;
  Matrix<Double> feed_offset(2,nRec);
  feed_offset = 0;
  Matrix<Complex> feed_response(nRec,nRec);
  feed_response = Complex(0.0,0.0);
  for (Int rec=0; rec<nRec; rec++) {
    feed_response(rec,rec) = Complex(1.0,0.0);
  }
  Vector<String> feed_type(nRec);
  feed_type(0) = "X";
  feed_type(1) = "Y";
  Vector<Double> feed_pos(3);
  feed_pos = 0.0;
  Vector<Double> feed_angle(nRec);
  feed_angle = -C::pi_4;                      // 0 for parallel dipoles
  // Fill the FEED subtable.
  MSFeed ms_feed = its_ms->feed();
  MSFeedColumns ms_feed_col(ms_feed);
  ms_feed.addRow (its_nantennas);
  for (Int i=0; i<its_nantennas; i++) {
    ms_feed_col.antennaId().put (i, i);
    ms_feed_col.feedId().put (i, 0);
    ms_feed_col.spectralWindowId().put (i, -1);
    //ms_feed_col.time().put (i, its_start_time + itsNrTimes*itsTimeStep/2.);
    ms_feed_col.time().put (i, its_start_time );
    //ms_feed_col.interval().put (i, itsNrTimes*itsTimeStep);
    ms_feed_col.interval().put (i, 0.);
    //std::cerr<<itsNrTimes*itsTimeStep<<std::endl;
    ms_feed_col.beamId().put (i, -1);
    ms_feed_col.beamOffset().put (i, feed_offset);
    ms_feed_col.polarizationType().put (i, feed_type);
    ms_feed_col.polResponse().put (i, feed_response);
    ms_feed_col.position().put (i, feed_pos);
    ms_feed_col.receptorAngle().put (i, feed_angle);
    ms_feed_col.numReceptors().put (i, 2);
  }
  ms_feed.flush();
}

void mscreate::fill_observation()
{
  MSObservation ms_obs = its_ms->observation();
  MSObservationColumns ms_obs_col(ms_obs);
  Vector<String> corrSchedule(1);
  corrSchedule = "corrSchedule";
  Vector<Double> timeRange(2);
  timeRange(0) = its_start_time;
  timeRange(1) = its_start_time ;
  // Data is public one year after end of observation.
  Double releaseDate = timeRange(1) ;
  // Fill the columns
  ms_obs.addRow();
  ms_obs_col.telescopeName().put (0, "21CMA");
  ms_obs_col.timeRange().put (0, timeRange);
  ms_obs_col.observer().put (0, "mscreate");
  ms_obs_col.scheduleType().put (0, "21CMA");
  ms_obs_col.schedule().put (0, corrSchedule);
  ms_obs_col.project().put (0, "mscreate");
  ms_obs_col.releaseDate().put (0, releaseDate);
  ms_obs_col.flagRow().put (0, False);
  ms_obs.flush();
}

void mscreate::fill_processor()
{
  MSProcessor ms_proc = its_ms->processor();
  MSProcessorColumns ms_proc_col(ms_proc);
  // Fill the columns
  ms_proc.addRow();
  ms_proc_col.type().put (0, "CORRELATOR");
  ms_proc_col.subType().put (0, "");
  ms_proc_col.typeId().put (0, -1);
  ms_proc_col.modeId().put (0, -1);
  ms_proc_col.flagRow().put (0, False);
  ms_proc.flush();
}

void mscreate::fill_state()
{
  MSState ms_state = its_ms->state();
  MSStateColumns ms_state_col(ms_state);
  // Fill the columns
  ms_state.addRow();
  ms_state_col.sig().put (0, True);
  ms_state_col.ref().put (0, False);
  ms_state_col.cal().put (0, 0.);
  ms_state_col.load().put (0, 0.);
  ms_state_col.subScan().put (0, 0);
  ms_state_col.obsMode().put (0, "");
  ms_state_col.flagRow().put (0, False);
  ms_state.flush();
}

void mscreate::update_times()
{
  // Calculate the interval, end, and central time.
  //Double interval = itsNrTimes*itsTimeStep;
  Double interval=its_end_time-its_start_time;
  //Double endTime = its_start_time + interval;
  Double midTime = (its_start_time + its_end_time) / 2;
  // Update all rows in FEED subtable.
  {
    MSFeed ms_sub = its_ms->feed();
    MSFeedColumns ms_sub_col(ms_sub);
    Vector<Double> val(ms_sub.nrow());
    val = midTime;
    ms_sub_col.time().putColumn (val);
    val = interval;
    ms_sub_col.interval().putColumn (val);
  }
  // Update all rows in POINTING subtable.
  {
    MSPointing ms_sub = its_ms->pointing();
    MSPointingColumns ms_sub_col(ms_sub);
    Vector<Double> val(ms_sub.nrow());
    val = midTime;
    ms_sub_col.time().putColumn (val);
    val = interval;
    ms_sub_col.interval().putColumn (val);
  }
  // Update all rows in OBSERVATION subtable.
  {
    MSObservation ms_obs = its_ms->observation();
    MSObservationColumns ms_obs_col(ms_obs);
    Vector<Double> timeRange(2);
    timeRange(0) = its_start_time;
    timeRange(1) = its_end_time;
    for (uInt i=0; i<ms_obs.nrow(); i++) {
      ms_obs_col.timeRange().put (i, timeRange);
    }
  }
}

void mscreate::write_time_step(raw_data_source& rds)//t in UTC in sec
{
  double t=rds.time();
  assert(t>its_end_time);
  
  Double time=(t+its_end_time)/2.0;
  Quantity qtime(time, "s");
  Double dt=t-its_end_time;
  its_end_time=t;
  cerr<<MEpoch(qtime, MEpoch::UTC)<<endl;
  its_frame->set (MEpoch(qtime, MEpoch::UTC));
  for (int field=0; field<its_nfields; ++field)
    {
      its_frame->set ((*its_phase_dir)[field]);
      vector<Vector<Double> > antuvw(its_nantennas);
      for (int j=0; j<its_nantennas; ++j) {
	MBaseline& mbl = (*its_ant_bl)[j];
	mbl.getRefPtr()->set(*its_frame);      // attach frame
	MBaseline::Convert mcvt(mbl, MBaseline::J2000);
	MVBaseline bas = mcvt().getValue();
	MVuvw jvguvw(bas, (*its_phase_dir)[field].getValue());
	antuvw[j] = Muvw(jvguvw, Muvw::J2000).getValue().getVector();
      }
      
      Vector<double> myuvw(3);
      for (int band=0; band<its_nbands; ++band)
	{
	  //defData = Complex(float(band), float(itsNrTimes));
	  IPosition shape(2, (*its_npol)[band], (*its_nchan)[band]);
	  //Array<Bool> defFlags(shape);
	  //defFlags = False;
	  //defData = Complex();
	  //int nrbasel = its_nantennas*(its_nantennas-1)/2;
	  //if (its_write_auto_corr) {
	  //  nrbasel += its_nantennas;
	  //}
	  int nrbasel=rds.num_of_baselines();
	  
	  // Add the number of rows needed.
	  int row_number = its_ms->nrow();
	  its_ms->addRow (nrbasel);
	  //std::cerr<<row_number<<" "<<its_ms->nrow()<<std::endl;
	  //Array<Float> sigma(IPosition(1, shape(0)));
	  //sigma = 1;
	  //Array<Float> weight(IPosition(1, shape(0)));
	  //weight = 1;
	  
	  
	  //defData = Complex(0,0);
	  //for (int j=0; j<its_nantennas; ++j)
	  //{
	  //int st = (its_write_auto_corr ? j : j+1);
	  //for (int i=st; i<its_nantennas; ++i)
	  for(int bl=0;bl<rds.num_of_baselines();++bl)
	    {
	      std::pair<int,int> antenna_pair(rds.antenna_pair(bl));
	      Array<Complex> defData(rds.data(field,band,bl));
	      Array<Bool> defFlags(rds.flags(field,band,bl));
	      Array<Float> sigma(rds.sigma(field,band,bl));
	      Array<Float> weight(pow(sigma,-2.0));
	      //uvw calculation
	      //see
	      //https://casa.nrao.edu/Memos/CoordConvention.pdf
	      //for uvw convension
	      
	      myuvw = antuvw[antenna_pair.first] - antuvw[antenna_pair.second];

#if 0
	      if(correct_w)
		{
		  double w=myuvw[2];
		  double v=myuvw[1];
		  double u=myuvw[0];
		  for(int i=0;i<ch_freq_vectors.at(band).size();++i)
		    {
		      constexpr double c=299792458.0;
		      constexpr double pi=std::atan(1)*4;
		      double freq=ch_freq_vectors.at(band).at(i);
		      double l=c/freq;
		      //std::cerr<<w/l<<endl;
		      IPosition p(2,0,i);
		      defData(p)*=std::exp(-std::complex<double>(0,1)*w/l*2.0*pi);
		      //defData(p)=std::exp(-std::complex<double>(0,w/l*2.0*pi));
		    }
		}
#endif
	      
	      
	      its_ms_col->data().put(row_number, defData);
	      its_ms_col->flag().put(row_number, defFlags);
	      its_ms_col->flagRow().put (row_number, False);
	      its_ms_col->time().put (row_number, time);
	      its_ms_col->antenna1().put (row_number, antenna_pair.first);
	      its_ms_col->antenna2().put (row_number, antenna_pair.second);
	      its_ms_col->feed1().put (row_number, 0);
	      its_ms_col->feed2().put (row_number, 0);
	      its_ms_col->dataDescId().put (row_number, band);
	      its_ms_col->processorId().put (row_number, 0);
	      its_ms_col->fieldId().put (row_number, field);
	      its_ms_col->interval().put (row_number, dt);
	      its_ms_col->exposure().put (row_number, dt);
	      its_ms_col->timeCentroid().put (row_number, time);
	      its_ms_col->scanNumber().put (row_number, 0);
	      its_ms_col->arrayId().put (row_number, 0);
	      its_ms_col->observationId().put (row_number, 0);
	      its_ms_col->stateId().put (row_number, 0);
	      its_ms_col->uvw().put (row_number, myuvw);
	      its_ms_col->weight().put (row_number, weight);
	      its_ms_col->sigma().put (row_number, sigma);
	      row_number++;
	    }
	}
    }
}

void mscreate::fill_baselines()
{
  MSAntenna ms_ant(its_ms->antenna());
  MSAntennaColumns ms_ant_col(ms_ant);
  Int nr = ms_ant.nrow();
  its_ant_bl = new Block<MBaseline>(nr);
  for (uInt j=0; j<ms_ant.nrow(); j++) {
    Vector<Double> antpos = ms_ant_col.position()(j);
    MVPosition blpos(antpos(0), antpos(1), antpos(2));
    (*its_ant_bl)[j] = MBaseline (MVBaseline(blpos), MBaseline::ITRF);
  }
}

void mscreate::flush()
{
  its_ms->flush(true);
}

casacore::Vector<casacore::Double> mscreate::calc_uvw(casacore::MBaseline bl,double utc_t,double ra,double dec)
{
  MPosition array_pos(casacore::MVPosition(0,0,0),MPosition::ITRF);
  MeasFrame frame(array_pos);
  Quantity qtime(utc_t,"s");
  frame.set(MEpoch(qtime,MEpoch::UTC));
  MDirection indir(MVDirection(Quantity(ra,"rad"),Quantity(dec,"rad")),MDirection::J2000);
  frame.set(indir);
  bl.getRefPtr()->set(frame);
  MBaseline::Convert mcvt(bl,MBaseline::J2000);
  MVBaseline bas=mcvt().getValue();
  MVuvw mvuvw(bas,indir.getValue());
  return Muvw(mvuvw,Muvw::J2000).getValue().getVector();
}
