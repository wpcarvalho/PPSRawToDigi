import pad_layout;
import root;

//string constraintsType = "homogeneous";
string constraintsType = "fixedDetectors";

string constraintNumber = "basic";
//string constraintNumber = "extended";

string analysisDir = "analysis/" + constraintsType + "." + constraintNumber;

string iteration = "5";

//---------------------------------------------------------------------------------------------------------------------

dotfactor = 10;

pen markupColor = paleyellow;

for (int i = 2; i <= 2; ++i) {
	ResetPads();
	
	if (constraintNumber != "basic") {
		NewPad(drawAxes = false);
		picture p; 
		label(p, rotate(0)*Label("V detectors"));
		attach(bbox(p, 1mm, nullpen, Fill(markupColor)));
		
		NewPad(drawAxes = false);
		picture p; 
		label(p, rotate(0)*Label("U detectors"));
		attach(bbox(p, 1mm, nullpen, Fill(markupColor)));
		NewRow();
	}

	pad pShRV = NewPad("detector number", "shift in $v\quad(\rm\mu m)$");
	pad pShRU = NewPad("detector number", "shift in $u\quad(\rm\mu m)$");
	if (constraintNumber != "basic")
		NewRow();
	pad pRotZV = NewPad("detector number", "rotation around $z\quad(\rm mrad)$");
	pad pRotZU;
	if (constraintNumber != "basic")
		pRotZU = NewPad("detector number", "rotation around $z\quad(\rm mrad)$");

	// ----- OPTICAL RESULTS ---

	string filename = "../" + analysisDir + "/RP" + ((string) i) + "/optical/results_Ideal.xml";
	write(filename);

	file f = input(filename);
	while (!eof(f)) {
		string line = f;
		string[] bits = split(line, "\"");
		//write(line);

		int id = -1;
		real sh_r = 0, sh_r_e = 0, rot_z = 0, rot_z_e = 0;

		for (int j = 0; j < bits.length; ++j) {
			//write("> ", bits[j]);
			if (find(bits[j], "id=") >= 0) id = (int) bits[++j];
			if (find(bits[j], "sh_r=") >= 0) sh_r = (real) bits[++j];
			if (find(bits[j], "sh_r_e=") >= 0) sh_r_e = (real) bits[++j];
			if (find(bits[j], "rot_z=") >= 0) rot_z = (real) bits[++j];
			if (find(bits[j], "rot_z_e=") >= 0) rot_z_e = (real) bits[++j];
		}
	
		if (id < 0)
			continue;
	
		int det = id % 10;

		sh_r = -sh_r;
		sh_r_e = 5.;

		rot_z_e = 0.5;

		SetPad((det % 2 == 0) ? pShRV : pShRU);
		draw((det, sh_r), marker(scale(2.5)*polygon(3), red, Fill));
		draw((det, sh_r)--(det, sh_r + sh_r_e), red, EndBar);
		draw((det, sh_r)--(det, sh_r - sh_r_e), red, EndBar);

		SetPad((det % 2 == 0 || constraintNumber == "basic") ? pRotZV : pRotZU);
		draw((det, rot_z), marker(scale(2.5)*polygon(3), red, Fill));
		draw((det, rot_z)--(det, rot_z + rot_z_e), red, EndBar);
		draw((det, rot_z)--(det, rot_z - rot_z_e), red, EndBar);
	}
	
	// ----- TRACK-BASED RESULTS ------

	string filename = "../" + analysisDir + "/RP" + ((string) i) + "/iteration" + iteration + "/cumulative_Jan.xml";
	write(filename);

	file f = input(filename);
	while (!eof(f)) {
		string line = f;
		string[] bits = split(line, "\"");
		//write(line);

		int id = -1;
		real sh_r = 0, sh_r_e = 0, rot_z = 0, rot_z_e = 0;

		for (int j = 0; j < bits.length; ++j) {
			//write("> ", bits[j]);
			if (find(bits[j], "id=") >= 0) id = (int) bits[++j];
			if (find(bits[j], "sh_r=") >= 0) sh_r = (real) bits[++j];
			if (find(bits[j], "sh_r_e=") >= 0) sh_r_e = (real) bits[++j];
			if (find(bits[j], "rot_z=") >= 0) rot_z = (real) bits[++j];
			if (find(bits[j], "rot_z_e=") >= 0) rot_z_e = (real) bits[++j];
		}
	
		//write("sh_r = ", sh_r);
	
		if (id < 0)
			continue;
	
		int det = id % 10;

		SetPad((det % 2 == 0) ? pShRV : pShRU);
		if (sh_r_e > 0) {
			draw((det, sh_r), marker(scale(2.5)*unitcircle, blue, Fill));
			draw((det, sh_r)--(det, sh_r + sh_r_e), blue, EndBar);
			draw((det, sh_r)--(det, sh_r - sh_r_e), blue, EndBar);
		} else {
			draw((det, sh_r), marker(scale(3.5)*polygon(4), Fill));
		}

		SetPad((det % 2 == 0 || constraintNumber == "basic") ? pRotZV : pRotZU);
		if (rot_z_e > 0) {
			draw((det, rot_z), marker(scale(2.5)*unitcircle, blue, Fill));
			draw((det, rot_z)--(det, rot_z + rot_z_e), blue, EndBar);
			draw((det, rot_z)--(det, rot_z - rot_z_e), blue, EndBar);
		} else {
			draw((det, sh_r), marker(scale(3.5)*polygon(4), Fill));
		}
	}


	SetPad(pShRV);
	limits((-0.5, -50), (9.5, +50), Crop);

	SetPad(pShRU);
	limits((-0.5, -50), (9.5, +50), Crop);

	SetPad(pRotZV);
	limits((-0.5, -10), (9.5, +10), Crop);

	SetPad(pRotZU);
	limits((-0.5, -10), (9.5, +10), Crop);
	
	NewRow();
	NewPad("angle   $(\rm mrad)$", "");
	string filename = "../" + analysisDir + "/RP" + ((string) i) + "/iteration" + iteration + "/statistics.root";
	draw(xscale(1e3), rGetObj(filename, "fitAxHist_selected"), "", heavygreen+1pt, "$a_x$");
	Legend bla;
	bla.operator init("entries = $" + format("%.0f", robj.rExec("GetEntries")) + "$", p=invisible);
	currentpicture.legend.push(bla);	
	Legend bla;
	bla.operator init("RMS = $" + format("%.1f", robj.rExec("GetRMS")*1e3) + "\ \rm mrad$", p=invisible);
	currentpicture.legend.push(bla);
	draw(xscale(1e3), rGetObj(filename, "fitAyHist_selected"), "", magenta+1pt, "$a_y$");
	Legend bla;
	bla.operator init("RMS = $" + format("%.1f", robj.rExec("GetRMS")*1e3) + "\ \rm mrad$", p=invisible);
	currentpicture.legend.push(bla);
	xlimits(-50, +50, Crop);
	frame legendFrame = legend();

	NewPad(drawAxes = false, autoSize=false);
	attach(legendFrame, (-100, 0));

	string outFilename = constraintsType + "." + constraintNumber + "/RP" + ((string) i);
	write("--> ", outFilename);
	GShipout(outFilename, NW);
}
