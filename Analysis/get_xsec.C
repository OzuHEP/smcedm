void get_xsec(TString infile_name, TString tree_name)
{
	auto file = new TFile(infile_name);
	TTreeReader myReader(tree_name, file);
	auto n_events = myReader.GetEntries(1);
	double xsec = 0.;
	//Get cross-section
	TTreeReaderArray<float>  ra_event_weight(myReader, "Weight.Weight");

	for (int i_event = 0; i_event < n_events; ++i_event){
		myReader.SetEntry(i_event);
		xsec += ra_event_weight.At(0);
		//std::cout << "i_event:" << i_event << "\tcross-section:" << ra_event_weight.At(0) << std::endl;
	}
	std::cout << "total cross-section:" << xsec << std::endl;
}