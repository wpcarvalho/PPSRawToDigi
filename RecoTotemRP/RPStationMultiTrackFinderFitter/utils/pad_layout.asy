/****************************************************************************
*
* This is a part of the "ROOT extension to Asymptote" project.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*
* \file pad_layout.asy
* \brief Macros for an intelligent layout of drawings into a final figure.
****************************************************************************/

import graph;
from settings access verbose;

//----------------------------------------------------------------------------------------------------
//------------------------------------- standard settings --------------------------------------------
//----------------------------------------------------------------------------------------------------

pen std_pens[] = { black, red, blue, heavygreen, cyan, magenta, orange, green, olive };
pen std_colors[] = std_pens;
pen stdPens[] = std_pens;

pen StdPen(int i, int base=0)
{
	if (base == 0)
		base = std_pens.length;

	int group = quotient(i, base);
	int idx = i % base;

	pen p = solid;
	if (group == 0) p = solid;
	if (group == 1) p = dashed;
	if (group == 2) p = dashdotted;
	if (group == 3) p = longdashed;

	return p + std_pens[idx];
}

pair O = (0, 0);


//----------------------------------------------------------------------------------------------------
//-------------------------------------- pad class and routines --------------------------------------
//----------------------------------------------------------------------------------------------------

/// default pad size
real xSizeDef = 6cm, ySizeDef = 6cm;

/// default axis ticks
ticks xTicksDef = LeftTicks, yTicksDef = RightTicks;

/// default grid hints
int xGridHintDef = 0, yGridHintDef = 0;

/**
 *\brief The class for a pad, a building block of a figure.
 *
 * New instace are to be created by one of the NewPad functions.
 **/
struct pad
{
	/// the containter for graphics in this pad
	picture pic;

	/// whether to add axes to this pad
	bool drawAxes = true;

	/// whether axes shall be placed above other drawings
	bool axesAbove = false;

	/// axes labels
	string xLabel = "", yLabel = "";

	/// axes tick styles
	ticks xTicks = xTicksDef, yTicks = yTicksDef;

	/// whether to apply automatic sizing
	bool autoSize = true;

	/// size of this pad
	real xSize = xSizeDef, ySize = ySizeDef;

	/// hint for placement within a grid of pads via GShipout
	int xGridHint, yGridHint;

	/// TODO
	bool fixed = false;

	/// TODO
	real x, y;
};

//----------------------------------------------------------------------------------------------------

/// the reference to the current/active pad
pad currentpad;

///\brief the collection of pads
/// empty in the begging => NewPad must be called in the begging
pad[] pad_collection = {};

//----------------------------------------------------------------------------------------------------

/**
 *\brief The common base for creating new pads. No need to use it directly (use NewPad functions instead). 
 **/
pad NewPadBase(	bool drawAxes,
				bool axesAbove,
				string xLabel, string yLabel,
				ticks xTicks, ticks yTicks,
				bool autoSize,
				real xSize, real ySize)
{
	if (verbose > 0) {
		write(">> NewPadBase");
		write("    drawAxes = ", drawAxes);
		write("    axesAbove = ", axesAbove);
		write("    xLabel = ", xLabel);
		write("    yLabel = ", yLabel);
		//write("    xTicks = ", xTicks);
		//write("    yTicks = ", yTicks);
		write("    autoSize = ", autoSize);
		write("    xSize = ", xSize);
		write("    ySize = ", ySize);
	}

	currentpad = new pad;
	currentpicture = currentpad.pic;
	pad_collection.push(currentpad);

	currentpad.drawAxes = drawAxes;
	currentpad.axesAbove = axesAbove;
	currentpad.xLabel = xLabel;
	currentpad.yLabel = yLabel;
	currentpad.xTicks = xTicks;
	currentpad.yTicks = yTicks;
	currentpad.autoSize = autoSize;
	currentpad.xSize = xSize;
	currentpad.ySize = ySize;
	
	return currentpad;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Creates new pad with automatic grid hints. 
 **/
pad NewPad(	bool drawAxes = true,
			bool axesAbove = false,
			string xLabel = "", string yLabel = "",
			ticks xTicks = xTicksDef, ticks yTicks = yTicksDef,
			bool autoSize = true,
			explicit real xSize = xSizeDef, explicit real ySize = ySizeDef)
{
	if (verbose > 0)
		write(">> NewPad, auto grid hints");

	NewPadBase(drawAxes, axesAbove, xLabel, yLabel, xTicks, yTicks, autoSize, xSize, ySize);
	
	currentpad.xGridHint = xGridHintDef;
	currentpad.yGridHint = yGridHintDef;
	++xGridHintDef;

	return currentpad;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Creates new pad with manual grid hints.
 **/
pad NewPad(	bool drawAxes = true,
			bool axesAbove = false,
			string xLabel = "", string yLabel = "",
			ticks xTicks = xTicksDef, ticks yTicks = yTicksDef,
			bool autoSize = true,
			explicit real xSize = xSizeDef, explicit real ySize = ySizeDef,
			explicit int xGridHint, int yGridHint)
{
	if (verbose > 0)
		write(">> NewPad, manual grid hints: ", xGridHint, yGridHint);

	NewPadBase(drawAxes, axesAbove, xLabel, yLabel, xTicks, yTicks, autoSize, xSize, ySize);
	
	currentpad.xGridHint = xGridHint;
	currentpad.yGridHint = yGridHint;
	
	return currentpad;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Sets x and y axes labels.
 **/
void SetPadLabels(pad p=currentpad, string xl, string yl)
{
	p.xLabel = xl;
	p.yLabel = yl;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Sets grid hint for the given pad.
 **/
void SetGridHint(pad p = currentpad, int x, int y)
{
	p.xGridHint = x;
	p.yGridHint = y;

	if (verbose > 0)
		write(">> SetGridHint: ", x, y);
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief TODO
 **/
void FixPad(pad p=currentpad, real x, real y)
{
	currentpad.fixed = true;
	currentpad.x = x;
	currentpad.y = y;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Starts new row.
 **/
void NewRow()
{
	++yGridHintDef;
	xGridHintDef = 0;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Starts new page.
 **/
void NewPage()
{
	if (verbose > 0)
		write(">> NewPage");

	newpage();
	xGridHintDef = yGridHintDef = 0;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Clears and removes all pads
 **/
void ResetPads()
{
	if (verbose > 0)
		write(">> ResetPads");

	pad_collection.delete();
	xGridHintDef = yGridHintDef = 0;
}

ResetPads();

//----------------------------------------------------------------------------------------------------

/**
 *\brief Sets current/active pad
 **/
void SetPad(pad p)
{
	currentpad = p;
	currentpicture = p.pic;
}

//----------------------------------------------------------------------------------------------------

void write(pad p)
{
	write("xLabel='"+p.xLabel+"', yLabel='"+p.yLabel+"', fixed=", p.fixed);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------- shipout routines -----------------------------------------
//----------------------------------------------------------------------------------------------------

/**
 *\brief Exports an assembled figure to a file. No need to call this function directly.
 **/
void FinalShipout(string prefix, picture final, real margin, pen p, filltype filltype)
{
	frame ff = bbox(final, margin, p, filltype);
	if (verbose > 0) {
		write(">> final size in cm:");
		write("  width = ", (max(ff).x - min(ff).x) / 1cm);
		write("  height = ", (max(ff).y - min(ff).y) / 1cm);
	}
	shipout(prefix, ff);

	ResetPads();
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Applies final modifications (adds axes, sizing, ...) to a pad before shipout. No need to call this function directly.
 **/
void FinalizePad(pad p)
{
	if (p.drawAxes)
	{
		xaxis(p.pic, Label(p.xLabel, 1), BottomTop, p.xTicks, above=p.axesAbove);
		yaxis(p.pic, Label(p.yLabel, 1), LeftRight, p.yTicks, above=p.axesAbove);
	}
	
	pair min = p.pic.userMin();
	pair max = p.pic.userMax();

	if (p.autoSize) {
		size(p.pic, p.xSize, p.ySize, min, max); 
	}
}

//----------------------------------------------------------------------------------------------------

frame PadToFrame(pad p, pair alignment)
{
	pair p_sh = point(p.pic, alignment);
	//dot(p.pic, p_sh, heavygreen+7pt);
	p_sh = p.pic.calculateTransform() * p_sh;
	return shift(-p_sh) * p.pic.fit();
	
	frame f = shift(-p_sh) * bbox(p.pic, Fill(paleyellow));
	return f;
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Aligns pads to a grid, using GridHints.
 *
 *\param alignment The position of the alignment control point.
 *\param hSkip The separation between columns.
 *\param vSkip The separation between rows.
 *\param margin The margin size around the final figure.
 *\param p The pen to draw a frame around the final figure.
 *\param filltype The bounding box around the final figure is filled with this type.
 **/
void GShipout(string prefix=defaultfilename, pair alignment=(0, 0), real hSkip=1cm, real vSkip=1cm, real margin=5mm, 
	pen p = nullpen, filltype filltype = Fill(white))
{
	if (verbose > 0)
		write(">> GShipout: prefix = ", prefix);

	if (verbose > 0)
		write(">> pads total: ", pad_collection.length);

	// create list of pads to be exported (with no gaps in indexes)
	pad[] ePads;
	for (pad p : pad_collection) {
		if (!p.pic.empty() && !p.fixed)
			ePads.push(p);
	}

	if (verbose > 0)
		write(">> pads to shipout: ", ePads.length);

	// determine grid boundaries
	int col_min = 1, col_max = 0, row_min = 1, row_max = 0;
	for (int i = 0; i < ePads.length; ++i) {
		if (ePads[i].xGridHint < col_min) col_min = ePads[i].xGridHint;
		if (ePads[i].xGridHint > col_max) col_max = ePads[i].xGridHint;
		if (ePads[i].yGridHint < row_min) row_min = ePads[i].yGridHint;
		if (ePads[i].yGridHint > row_max) row_max = ePads[i].yGridHint;
	}
	int cols = col_max - col_min + 1;
	int rows = row_max - row_min + 1;

	if (cols < 1 || rows < 1) {
		write("WARNING in GShipout > Nothing to shipout");	
		return;
	}

	if (verbose > 0) {
		write(">> grid boundaries");
		write("  col_min = ", col_min);
		write("  col_max = ", col_max);
		write("  row_min = ", row_min);
		write("  row_max = ", row_max);
	}

	// determine row and column sizes
	real[] wR, wL, hT, hB;
	frame [] frames;
	for (int i = 0; i < ePads.length; ++i) {
		FinalizePad(ePads[i]);
		frame f = PadToFrame(ePads[i], alignment);
		frames.push(f);

		pair max = max(f), min = min(f);
		int col = ePads[i].xGridHint - col_min, row = ePads[i].yGridHint - row_min;
		if (verbose > 0) {
			write("* frame ", i);
			write("  cell: ", row, col);
			write("  min: ", min);
			write("  max: ", max);
		}

		if (! wR.initialized(col)) {
			wR[col] = max.x;
			wL[col] = -min.x;
		} else {
			if (wR[col] < max.x) wR[col] = max.x;
			if (wL[col] < -min.x) wL[col] = -min.x;
		}

		if (! hT.initialized(row)) {
			hT[row] = max.y;
			hB[row] = -min.y;
		} else {
			if (hT[row] < max.y) hT[row] = max.y;
			if (hB[row] < -min.y) hB[row] = -min.y;
		}
	}

	if (verbose > 0) {
		write(">> column widths: right, left");
		for (int i : wR.keys)
			write(i, wR[i], wL[i]);
		write(">> row heights: top, bottom");
		for (int i : hT.keys)
			write(i, hT[i], hB[i]);
	}

	// calculate offsets
	real[] xOffsets, yOffsets;
	real x = 0, y = 0;
	if (verbose > 0)
		write(">>  x offsets");
	for (int c = 0; c < cols; ++c) {
		if (c > 0) {
			if (wR.initialized(c - 1)) x += wR[c - 1];
			x += hSkip;
			if (wL.initialized(c)) x += wL[c];
		}
		xOffsets[c] = x;
		if (verbose > 0)
			write("  ", c, x);
	}

	if (verbose > 0)
		write(">>  y offsets");
	for (int r = 0; r < rows; ++r) {
		if (r > 0) {
			if (hB.initialized(r - 1)) y -= hB[r - 1];
			y -= vSkip;
			if (hT.initialized(r)) y -= hT[r];
		}
		yOffsets[r] = y;
		if (verbose > 0)
			write("  ", r, y);
	}

	// align the non-fixed pads
	if (verbose > 0)
		write(">>  final alignment");
	picture final;
	for (int i = 0; i < frames.length; ++i) {
		pair pos = (xOffsets[ePads[i].xGridHint - col_min], yOffsets[ePads[i].yGridHint - row_min]);
		if (verbose > 0)
			write("  pad " + format("%u", i) + ": ", pos);
		add(final, frames[i], position=pos);
	}
	
	// add fixed pads
	for (pad p : pad_collection) {
		if (!p.pic.empty() && p.fixed) {
			FinalizePad(p);
			add(final, PadToFrame(p, alignment), position=(p.x, p.y));
		}
	}

	// finish
	FinalShipout(prefix, final, margin, p, filltype);
}

//----------------------------------------------------------------------------------------------------

///\brief Alternative name to GShipout
void Shipout(string prefix=defaultfilename, pair alignment=(0, 0), real hSkip=1cm, real vSkip=1cm, real margin=1mm, 
	pen p = nullpen, filltype filltype = Fill(white)) = GShipout;


//----------------------------------------------------------------------------------------------------

/**
 *\brief The function that is automatically called at the exit time.
 **/
void exitfunction()
{
	if (needshipout())
		GShipout();
}

atexit(exitfunction);

//----------------------------------------------------------------------------------------------------
//-------------------------------------------- markers -----------------------------------------------
//----------------------------------------------------------------------------------------------------

struct mark {
	path shape = nullpath;
	pen color = black;
	real scale = 1.;
	bool fill = true;
}

mark mark(path sh = nullpath, pen co=black, real sc=1., bool fi=true)
{
	mark m;
	m.shape = sh;
	m.color = co;
	m.scale = sc;
	m.fill = fi;
	return m;
}

mark operator+ (mark m, path s)
{
	return mark(s, m.color, m.scale, m.fill);
}

mark operator+ (mark m, pen c)
{
	return mark(m.shape, c, m.scale, m.fill);
}

mark operator+ (mark m, real s)
{
	return mark(m.shape, m.color, s, m.fill);
}

mark operator+ (mark m, bool f)
{
	return mark(m.shape, m.color, m.scale, f);
}

mark operator* (transform t, mark m)
{
	return mark(t * m.shape, m.color, m.scale, m.fill);
}

void write(mark m)
{
	write("shape = ", m.shape);
	write("color = ", m.color);
	write("scale = ", m.scale);
	write("fill = ", m.fill);
}

marker operator cast(mark m)
{
	picture p;
	if (m.fill)
		filldraw(p, scale(m.scale)*m.shape, m.color, m.color);
	else
		filldraw(p, scale(m.scale)*m.shape, nullpen, m.color);
	return marker(p.fit());
}

path make_star(int nodes, real r_out=1., real r_in=0.4)
{
	path p;
	pair r = (0, 1);
	real de_a = 360./nodes/2;

	for (int i = 0; i < nodes; ++i) {
		p = p--r_out*r;
		r = rotate(de_a)*r;
		p = p--r_in*r;
		r = rotate(de_a)*r;
	}
	p = p--cycle;

	return p;
}

// predefined markers
real sin_60 = sqrt(3)/2;
mark mCi = mark(unitcircle);
mark mSq = mark(polygon(4));
mark mTU = mark(polygon(3));
mark mTL = mark(rotate(90)*mTU.shape);
mark mTD = mark(rotate(180)*mTU.shape);
mark mTR = mark(rotate(270)*mTU.shape);

mark mPl = mark((0,0)--(0,-1)--(0,0)--(1,0)--(0,0)--(0,1)--(0,0)--(-1,0)--cycle);
mark mCr = mark(rotate(45)*mPl.shape);

mark mPl2 = mark(scale(1./sqrt(10))*((1,1)--(1,3)--(-1,3)--(-1,1)--(-3,1)--(-3,-1)--
	(-1,-1)--(-1,-3)--(1,-3)--(1,-1)--(3,-1)--(3,1)--cycle));
mark mCr2 = mark(rotate(45)*mPl2.shape);

mark mSt4 = mark(make_star(4));
mark mSt5 = mark(make_star(5));
mark mSt6 = mark(make_star(6));
mark mSt = mSt5;

// standard mark collection
mark std_marks[] = { mCi, mSq, mTU, mTD, mSt };

mark StdMark(int i, int base=0)
{
	if (base == 0)
		base = std_marks.length;
	return std_marks[i % base];
}

//----------------------------------------------------------------------------------------------------
//--------------------------------------- legend routines --------------------------------------------
//----------------------------------------------------------------------------------------------------

/**
 *\brief Adds an item to the legend.
 **/
void AddToLegend(string label, pen p = invisible, frame f)
{
	Legend bla;
	bla.operator init(label, p, f);
	bla.plabel = black;
	currentpicture.legend.push(bla);
}

//----------------------------------------------------------------------------------------------------

void AddToLegend(string label, pen p = invisible, marker m = nomarker)
{
	AddToLegend(label, p, m.f);
}

//----------------------------------------------------------------------------------------------------

/**
 *\brief Alternative name to AddToLegend.
 **/
void AddLegend(string label, pen p = invisible, marker m = nomarker) = AddToLegend;

//----------------------------------------------------------------------------------------------------

/**
 *\brief Alternative name to AddToLegend.
 **/
void LegendAdd(string label, pen p = invisible, marker m = nomarker) = AddToLegend;

//----------------------------------------------------------------------------------------------------

/**
 *\brief Clears the legend.
 **/
void ClearLegend()
{
	currentpicture.legend.delete();
}

//----------------------------------------------------------------------------------------------------

pen legendLabelPen = black;

/**
 *\brief Converts a legend item into a picture.
 *\param lineLength The length of a line in the legend.
 *\param width If zero, the label will be formatted in one row. If non-zero, the label will be
 *             formatted in a block of that width.
 **/
picture LegendItem(Legend item, real lineLength, real width=0)
{
	picture pic;
	pair z1 = (0, 0);
	real gap = 1mm;
	pair z2 = z1 + (lineLength+gap, 0);

	if (!item.above && !empty(item.mark))
		marknodes(pic, item.mark, interp(z1, z2, 0.5));

	if (lineLength > 0)
		draw(pic, z1--z2, item.p);
	
	if (item.above && !empty(item.mark))
		marknodes(pic, item.mark, interp(z1, z2, 0.5));

	string label = "\strut "+item.label+"\strut";
	if (width > 0) {
		label = "\vbox{\hsize"+format("%f", width)+"pt\noindent" + label + "}";
	}
	//label(pic, label, z2, E, item.plabel);
	label(pic, label, z2, E, legendLabelPen);
	return pic;
}

//----------------------------------------------------------------------------------------------------

/*
 *\brief Puts legend items into a column layout.
 */
picture LegendColumns(Legend[] items, int columns, real lineLength, real maxColWidth, real hSkip, real vSkip,
			bool stretch)
{
	int N = items.length;

	// make the columns number sane
	if (columns < 1)
		columns = 1;
	if (columns > N)
		columns = N;

	// adjust line length (in case of large Markers)
	/*
  	real maxSymbSize = 0;
	for (int i = 0; i < N; ++i) {
		real size = size(items[i].mark).x;
		if (items[i].p != invisible)
			size = max(size, lineLength);
		if (maxSymbSize < size)
			maxSymbSize = size;
	}
	lineLength = maxSymbSize;
	*/

	// make legend item pictures
	picture itemPic[];
	real incH[];
	incH[0] = 0.;
	for (int i = 0; i < N; ++i) {
		picture p = LegendItem(items[i], lineLength, maxColWidth);
		itemPic.push(p);
		incH[i+1] = incH[i] + size(p).y;
	}

	// find optimal column breaks
	int breaks[];
	real minH = +inf;
	if (columns > 1) {
		int B = columns - 1;

		//write("N = ", N);
		//write("columns = ", columns);
		//write("B = ", B);

		int b[];
		for (int i = 0; i < B; ++i)
			b.push(i+1);
		--b[B-1];
			
		while (true) {
			//write("--------------------");
			int l = B - 1;
			while (l >= 0) {
				if (b[l] <= N-(B-l))
					break;
				--l;
			}

			//write("l = ", l);

			if (l < 0)
				break;

			++b[l];
			for (int i = l+1; i < B; ++i)
				b[i] = b[i-1]+1;

			//write(b);

			real maxH = -inf;
			for (int i = 0; i < columns; ++i) {
				int b1 = (i > 0) ? b[i-1] : 0;
				int b2 = (i < B) ? b[i] : N;
				real h = incH[b2] - incH[b1] + (b2-b1)*vSkip;
				maxH = max(maxH, h);
			}
			
			if (maxH < minH) {
				//write("* maxH=", maxH);
				minH = maxH;
				breaks = copy(b);
				//write(breaks);
			}
		}
	} else {
		minH = incH[N];
	}

	//write("====================");
	//write("minH = ", minH);
	//write(breaks);

	// find stretches
	real stretches[];
	for (int i = 0; i < columns; ++i) {
		int b1 = (i > 0) ? breaks[i-1] : 0;
		int b2 = (i < columns-1) ? breaks[i] : N;
		real h = incH[b2] - incH[b1] + (b2-b1)*vSkip;

		stretches[i] = (b2-b1 > 1) ? (minH - h) / (b2-b1-1) : 0.;
	}
	
	//write("stretches = ", stretches);
	//write("====================");

	// align the legend pictures
	picture p;
	int c = 0;
	real maxW = 0, x = 0, y = 0;	
	for (int i = 0; i < N; ++i) {
		pair size = size(itemPic[i]);
		maxW = max(size.x, maxW);

		y -= size.y/2;
		add(p, itemPic[i], (x, y));
		//add(p, bbox(itemPic[i], nullpen, Fill(palered)), (x, y));
		real skip = vSkip;
		if (c < columns && stretch)
			skip += stretches[c];
		y -= size.y/2 + skip;	
		
		if (c < breaks.length && i == breaks[c] - 1) {
			y = 0;
			x += maxW + hSkip;
			maxW = 0;
			++c;
		}
	}

	return p;
}

//----------------------------------------------------------------------------------------------------

/*
 *\brief Puts the legend picture into a correctly aligned frame.
 */
frame BuildLegend(picture pic=currentpicture, string title="", int columns=1, real lineLength=1cm, real colWidth=0.,
		real hSkip=1mm, real vSkip=0mm, bool stretch=true, real xmargin=1mm, real ymargin=1mm,
		pen framePen=black, pair alignment=NE)
{
	//write("BuildLegend");
	//write("  title = ", title);

	picture lp = LegendColumns(pic.legend, columns, lineLength, colWidth, hSkip, vSkip, stretch);

	if (length(title) > 0)
		label(lp, "\strut "+title+"\strut", (0, vSkip), NE, legendLabelPen);
	
	// draw bounding box
	pair max = max(lp, true);
	pair min = min(lp, true);
	draw(lp, (min.x - xmargin, min.y - ymargin)--(max.x + xmargin, min.y - ymargin)--
	  (max.x + xmargin, max.y + ymargin)--(min.x - xmargin, max.y + ymargin)--cycle, framePen);
	
	// fit to frame such that the selected point is at (0, 0) in PostScript coordinates
	pair p = lp.calculateTransform() * truepoint(lp, alignment);
	frame F = shift(-p) * lp.fit();
	//frame F = shift(-p) * bbox(lp, 1mm, nullpen, Fill(paleblue));
	return F;
}

//----------------------------------------------------------------------------------------------------

/*
 *\brief The standard way to attach the legend to the current pad.
 * The legend is aligned such that legAlig point of the legend is placed
 * on top of the picAlig point of the picture.
*/
void AttachLegend(string title="", int columns=1, real colWidth=0., bool stretch=true, 
	pair legAlig=NE, pair picAlig=NE, filltype fillT = Fill(white))
{
	add(BuildLegend(title, columns, colWidth=colWidth, stretch=stretch, legAlig), point(picAlig), fillT);
}
