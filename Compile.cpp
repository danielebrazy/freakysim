void compile(const char option) {
	gROOT->LoadMacro("ParticleType.cpp+");
	gROOT->LoadMacro("ResonanceType.cpp+");
	gROOT->LoadMacro("Particle.cpp+");
	if (option == 's') gROOT->LoadMacro("Simulate.cpp+");
	else if (option == 't') gROOT->LoadMacro("test.cpp+");
	else if (option == 'a') gROOT->LoadMacro("Analyze.cpp+");
	else std::cout << "Invalid option.";
}

