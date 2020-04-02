/*
 * Visualization.cpp
 *
 *  Created on: Apr 28, 2012
 *      Author: mixlmay
 */

#include "Visualization.h"

inline double AbsoluteValue(double val) {
	if (val < 0.) return(val*-1);
	else return(val);
}

inline double Sign(double val) {
	if (val < 0.) return (-1.);
	else return(1.);
}

Fl_Text_Display::Style_Table_Entry styletable[] = {     // Style table
  { FL_BLACK,      FL_COURIER,        FL_NORMAL_SIZE }, // A - Plain
  { FL_DARK_GREEN, FL_COURIER_ITALIC, FL_NORMAL_SIZE }, // B - Line comments
  { FL_DARK_GREEN, FL_COURIER_ITALIC, FL_NORMAL_SIZE }, // C - Block comments
  { FL_BLUE,       FL_COURIER,        FL_NORMAL_SIZE }, // D - Strings
  { FL_DARK_RED,   FL_COURIER,        FL_NORMAL_SIZE }, // E - Directives
  { FL_DARK_RED,   FL_COURIER_BOLD,   FL_NORMAL_SIZE }, // F - Types
  { FL_BLUE,       FL_COURIER_BOLD,   FL_NORMAL_SIZE }  // G - Keywords
};

// copied and modified from http://www.fltk.org/doc-1.1/editor.html
// 'style_parse()' - Parse text and produce style data.
//
const char         *code_keywords[] = { // List of known C/C++ keywords...
                   };
const char         *code_types[] = {    // List of known C/C++ types...
                   };


//
// 'compare_keywords()' - Compare two keywords...
//

extern "C" {
  int
  compare_keywords(const void *a,
                   const void *b) {
    return (strcmp(*((const char **)a), *((const char **)b)));
  }
}

void
style_parse(const char *text,
            char       *style,
            int        length) {
  char             current;
  int             col;
  int             last;
  char             buf[255],
             *bufptr;
  const char *temp;

  for (current = *style, col = 0, last = 0; length > 0; length --, text ++) {
    if (current == 'A') {
      // Check for directives, comments, strings, and keywords...
      if (col == 0 && *text == '#') {
        // Set style to directive
        current = 'E';
      } else if (strncmp(text, "1", 1) == 0) {
        current = 'B';
      } else if (strncmp(text, "2", 1) == 0) {
        current = 'E';
      } else if (strncmp(text, "\\\"", 2) == 0) {
        // Quoted quote...
        *style++ = current;
        *style++ = current;
        text ++;
        length --;
        col += 2;
        continue;
      } else if (*text == '\"') {
        current = 'D';
      } else if (!last && islower(*text)) {
        // Might be a keyword...
        for (temp = text, bufptr = buf;
             islower(*temp) && bufptr < (buf + sizeof(buf) - 1);
             *bufptr++ = *temp++);

        if (!islower(*temp)) {
          *bufptr = '\0';

          bufptr = buf;

          if (bsearch(&bufptr, code_types,
                      sizeof(code_types) / sizeof(code_types[0]),
                      sizeof(code_types[0]), compare_keywords)) {
            while (text < temp) {
              *style++ = 'F';
              text ++;
              length --;
              col ++;
            }

            text --;
            length ++;
            last = 1;
            continue;
          } else if (bsearch(&bufptr, code_keywords,
                             sizeof(code_keywords) / sizeof(code_keywords[0]),
                             sizeof(code_keywords[0]), compare_keywords)) {
            while (text < temp) {
              *style++ = 'G';
              text ++;
              length --;
              col ++;
            }

            text --;
            length ++;
            last = 1;
            continue;
          }
        }
      }
    } else if (current == 'C' && strncmp(text, "*/", 2) == 0) {
      // Close a C comment...
      *style++ = current;
      *style++ = current;
      text ++;
      length --;
      current = 'A';
      col += 2;
      continue;
    } else if (current == 'D') {
      // Continuing in string...
      if (strncmp(text, "\\\"", 2) == 0) {
        // Quoted end quote...
        *style++ = current;
        *style++ = current;
        text ++;
        length --;
        col += 2;
        continue;
      } else if (*text == '\"') {
        // End quote...
        *style++ = current;
        col ++;
        current = 'A';
        continue;
      }
    }

    // Copy style info...
    if (current == 'A' && (*text == '{' || *text == '}')) *style++ = 'G';
    else *style++ = current;
    col ++;

    last = isalnum(*text) || *text == '.';

    if (*text == '\n') {
      // Reset column and possibly reset the style
      col = 0;
      if (current == 'B' || current == 'E') current = 'A';
    }
  }
}


//
// 'style_init()' - Initialize the style buffer...
//

void
style_init(void) {
  char *style = new char[textbuf->length() + 1];
  const char *text = textbuf->text();

  memset(style, 'A', textbuf->length());
  style[textbuf->length()] = '\0';

  if (!stylebuf) stylebuf = new Fl_Text_Buffer(textbuf->length());

  style_parse(text, style, textbuf->length());

  stylebuf->text(style);
  delete[] style;
}


// copied and modified from http://www.fltk.org/doc-1.1/editor.html
// 'style_update()' - Update the style buffer...
//
void
style_update(int        pos,          // I - Position of update
             int        nInserted,    // I - Number of inserted chars
             int        nDeleted,     // I - Number of deleted chars
             int        nRestyled,    // I - Number of restyled chars
             const char *deletedText, // I - Text that was deleted
             void       *cbArg) {     // I - Callback data
  int  start,                         // Start of text
       end;                           // End of text
  char last,                          // Last style on line
       *style,                        // Style data
       *text;                         // Text data


  // If this is just a selection change, just unselect the style buffer...
  if (nInserted == 0 && nDeleted == 0) {
    stylebuf->unselect();
    return;
  }

  // Track changes in the text buffer...
  if (nInserted > 0) {
    // Insert characters into the style buffer...
    style = new char[nInserted + 1];
    memset(style, 'A', nInserted);
    style[nInserted] = '\0';

    stylebuf->replace(pos, pos + nDeleted, style);
    delete[] style;
  } else {
    // Just delete characters in the style buffer...
    stylebuf->remove(pos, pos + nDeleted);
  }

  // Select the area that was just updated to avoid unnecessary
  // callbacks...
  stylebuf->select(pos, pos + nInserted - nDeleted);

  // Re-parse the changed region; we do this by parsing from the
  // beginning of the line of the changed region to the end of
  // the line of the changed region...  Then we check the last
  // style character and keep updating if we have a multi-line
  // comment character...
  start = textbuf->line_start(pos);
  end   = textbuf->line_end(pos + nInserted - nDeleted);
  text  = textbuf->text_range(start, end);
  style = stylebuf->text_range(start, end);
  last  = style[end - start - 1];

  style_parse(text, style, end - start);

  stylebuf->replace(start, end, style);
  ((Fl_Text_Display*)cbArg)->redisplay_range(start, end);

  if (last != style[end - start - 1]) {
    // The last character on the line changed styles, so reparse the
    // remainder of the buffer...
    free(text);
    free(style);

    end   = textbuf->length();
    text  = textbuf->text_range(start, end);
    style = stylebuf->text_range(start, end);

    style_parse(text, style, end - start);

    stylebuf->replace(start, end, style);
    ((Fl_Text_Display *)cbArg)->redisplay_range(start, end);
  }

  free(text);
  free(style);
}

void style_unfinished_cb(int, void*) {
}


inline std::string IntToString(const unsigned int v){
	std::ostringstream output;
	if(!(output << v))
		throw BadConversion("Double2String could not convert value!");
	return(output.str());
}

unsigned int splitString(const std::string& str, std::vector<double>& data) {
	// will go into 'tokenize()' and 'removeComments()'
	unsigned int count_;
	unsigned int begin_;
	unsigned int end_;
	std::string delim_;

	double val;
	delim_ = " ,\t\f\n\v";

	//Point to the first token
	begin_ = str.find_first_not_of(delim_);
	end_ = str.find_first_of(delim_, begin_);
	count_ = 0;

	while (begin_ != std::string::npos){
		if (end_ != std::string::npos) {
			//data.push_back(str.substr(begin_, end_-begin_));
			val = boost::lexical_cast<double>(str.substr(begin_, end_-begin_));
			data.push_back(val);
			begin_ = str.find_first_not_of(delim_, end_);
			end_ = str.find_first_of(delim_, begin_);
			count_++;
		}
		else if (begin_ != std::string::npos && end_ == std::string::npos){
			//data.push_back(str.substr(begin_, str.length()-begin_));
			val = boost::lexical_cast<double>(str.substr(begin_, str.length()-begin_));
			data.push_back(val);
			begin_ = str.find_first_not_of(delim_, end_);
			count_++;
		}
	}

	return (count_);
}

void PlotBox::draw() {

	std::string buf;
	double TRANSLATE_X = (double)w()/2. - 300;
	double TRANSLATE_Y = (double)h()/2. + 100.;

	double SCALE_X = 30.;
	double SCALE_Y = -30.;

	//fl_clip(x(), y(), w(), h());
	//fl_push_matrix();
	// shift the origin to the lower left corner of the rectangle
	//fl_translate(w()/2, h()/2+100);
	// scale the AE data so that the axis span 20sec and 100dB
	//fl_scale(30,-30);
	fl_color(FL_DARK3);
	fl_rectf(x(), y(), w(), h());
	fl_font(FL_HELVETICA,15);
	for (unsigned int q = 1; q < ParticleData[step]->size();){
		fl_color(FL_RED);
		fl_circle((int)(TRANSLATE_X + SCALE_X * ParticleData[step]->at(q)),
				  (int)(TRANSLATE_Y + SCALE_Y * ParticleData[step]->at(q+1)),
				  (int)(              SCALE_X * ParticleData[step]->at(q+3)));
		fl_color(FL_DARK_RED);
		fl_pie((int)(TRANSLATE_X + SCALE_X * ParticleData[step]->at(q) - SCALE_X * ParticleData[step]->at(q+3)),
				  (int)(TRANSLATE_Y + SCALE_Y * ParticleData[step]->at(q+1) + SCALE_Y * ParticleData[step]->at(q+3)),
				  (int)(              SCALE_X * ParticleData[step]->at(q+3) * 2),
				  (int)(              SCALE_X * ParticleData[step]->at(q+3) * 2),
				  (int)(ParticleData[step]->at(q+2) * 57.3 + 45),
				  (int)(ParticleData[step]->at(q+2) * 57.3 + 135));
		fl_pie((int)(TRANSLATE_X + SCALE_X * ParticleData[step]->at(q) - SCALE_X * ParticleData[step]->at(q+3)),
				  (int)(TRANSLATE_Y + SCALE_Y * ParticleData[step]->at(q+1) + SCALE_Y * ParticleData[step]->at(q+3)),
				  (int)(              SCALE_X * ParticleData[step]->at(q+3) * 2),
				  (int)(              SCALE_X * ParticleData[step]->at(q+3) * 2),
				  (int)(ParticleData[step]->at(q+2) * 57.3 + 225),
				  (int)(ParticleData[step]->at(q+2) * 57.3 + 315));
		fl_color(FL_WHITE);
		buf = IntToString(q/8);
		fl_draw(buf.c_str(),
				(int)(TRANSLATE_X + SCALE_X * ParticleData[step]->at(q)),
				(int)(TRANSLATE_Y + SCALE_Y * ParticleData[step]->at(q+1)));
		q=q+8;
	}
	buf = "time = " + IntToString(step);
	fl_draw(buf.c_str(),
			x() + w()-100,
			y() + h()-10);
	for (unsigned int q = 1; q < ContactData[step]->size();){
		fl_color(FL_WHITE);
		if (ContactData[step]->at(q+8) < 0.0000000001 && ContactData[step]->at(q+9) < 0.0000000001) {
			fl_line_style(FL_DASH,1);
		}
		else {
			fl_line_style(FL_SOLID,
					(int) sqrt(
							ContactData[step]->at(q+8)*(int)ContactData[step]->at(q+8)
							+ ContactData[step]->at(q+9)*(int)ContactData[step]->at(q+9)) / 5);
		fl_line( (int)(TRANSLATE_X + SCALE_X * ParticleData[step]->at( 1 + (int)ContactData[step]->at(q)*8)),
				 (int)(TRANSLATE_Y + SCALE_Y * ParticleData[step]->at( 2 + (int)ContactData[step]->at(q)*8)),
				 (int)(TRANSLATE_X + SCALE_X * ParticleData[step]->at( 1 + (int)ContactData[step]->at(q+1)*8)),
				 (int)(TRANSLATE_Y + SCALE_Y * ParticleData[step]->at( 2 + (int)ContactData[step]->at(q+1)*8)));
		}
		q=q+10;
		fl_line_style(FL_SOLID,1);
	}
	//fl_pop_matrix();
	//fl_pop_clip();
}

void MyWindow::update_contact_display() {

	std::ostringstream output;
	double f_n, f_t, f, v_n, v_t, v;
	int p_1, p_2;

	for (unsigned int q = 1; q < this->MyPlotBox->ContactData[this->MyPlotBox->step]->size();){

		output.precision(4);

		p_1 = (int) this->MyPlotBox->ContactData[this->MyPlotBox->step]->at(q  );
		p_2 = (int) this->MyPlotBox->ContactData[this->MyPlotBox->step]->at(q+1);

		f_n = this->MyPlotBox->ContactData[this->MyPlotBox->step]->at(q+8);
		f_t = this->MyPlotBox->ContactData[this->MyPlotBox->step]->at(q+9);
		f   = sqrt(f_n*f_n + f_t*f_t);

		v_n = this->MyPlotBox->ContactData[this->MyPlotBox->step]->at(q+6);
		v_t = this->MyPlotBox->ContactData[this->MyPlotBox->step]->at(q+7);
		v   = sqrt(v_n*v_n * v_t*v_t)*Sign(v_n);

		if (f > 0.){
			if (v_t < -0.001 || v_t > 0.001)
				output << "2\t" << p_1 << "\t" << p_2 << "\t" << std::setw(10) << std::fixed << v << "\t" << f << std::endl;
			else
				output << "1\t" << p_1 << "\t" << p_2 << "\t" << std::setw(10) << std::fixed << v << "\t" << f << std::endl;
			textbuf->text(output.str().c_str());
			textbuf->append("\n");}
		q=q+10;
	}
}

inline void MyWindow::cb_the_one_backward_button_i() {
	if(this->MyPlotBox->step > 1) {
		this->MyPlotBox->step--;
		this->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		this->update_contact_display();
		}
}

inline void MyWindow::cb_the_one_forward_button_i() {
	if(this->MyPlotBox->step < this->MyPlotBox->max_step-2) {
		this->MyPlotBox->step++;
		this->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		this->update_contact_display();
		}
}

inline void MyWindow::cb_the_ten_backward_button_i() {
	if(this->MyPlotBox->step > 11) {
		this->MyPlotBox->step = this->MyPlotBox->step-10;
		this->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		this->update_contact_display();
		}
}

inline void MyWindow::cb_the_ten_forward_button_i() {
	if(this->MyPlotBox->step < this->MyPlotBox->max_step-12) {
		this->MyPlotBox->step = this->MyPlotBox->step+10;
		this->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		this->update_contact_display();
		}
}

inline void MyWindow::cb_the_hdr_backward_button_i() {
	if(this->MyPlotBox->step > 101) {
		this->MyPlotBox->step = this->MyPlotBox->step-100;
		this->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		this->update_contact_display();
		}
}

inline void MyWindow::cb_the_hdr_forward_button_i() {
	if(this->MyPlotBox->step < this->MyPlotBox->max_step-102) {
		this->MyPlotBox->step = this->MyPlotBox->step+100;
		this->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		this->update_contact_display();
		}
}

void cb_auto_update(void *v) {

	if(((MyWindow*)v)->MyPlotBox->step < ((MyWindow*)v)->MyPlotBox->max_step-2) {
		((MyWindow*)v)->MyPlotBox->step++;
		((MyWindow*)v)->MyPlotBox->redraw();
//			this->TheAeEventNrIndicator->value((double)tmp);
		((MyWindow*)v)->update_contact_display();
		}
	else {((MyWindow*)v)->MyPlotBox->step = 1; }

	Fl::repeat_timeout(0.4, cb_auto_update, v);
}

void MyWindow::cb_the_auto_forward_button(Fl_Widget *w, void *v) {

	if(((Fl_Light_Button*)w)->value() == 1)
		Fl::add_timeout(0.4,cb_auto_update, v);
	else
		Fl::remove_timeout(cb_auto_update);
}

MyWindow::MyWindow(int w, int h, const char* txt, std::string name_b) : Fl_Window(w,h,txt) {

	MyPlotBox = new PlotBox(10, 10, this->w()-320, this->h()-60, "");
	TheOneBackwardButton = new Fl_Button(this->w()/3-5-20,this->h()-35,30,30, "@<-");
	TheOneBackwardButton->callback(cb_the_one_backward_button, this);
	TheOneForwardButton = new Fl_Button(this->w()/3+5+20,this->h()-35,30,30, "@->");
	TheOneForwardButton->callback(cb_the_one_forward_button, this);
	TheTenBackwardButton = new Fl_Button(this->w()/3-10-80,this->h()-35,30,30, "@<-");
	TheTenBackwardButton->callback(cb_the_ten_backward_button, this);
	TheTenForwardButton = new Fl_Button(this->w()/3+10+80,this->h()-35,30,30, "@->");
	TheTenForwardButton->callback(cb_the_ten_forward_button, this);
	TheHdrBackwardButton = new Fl_Button(this->w()/3-15-120,this->h()-35,30,30, "@<-");
	TheHdrBackwardButton->callback(cb_the_hdr_backward_button, this);
	TheHdrForwardButton = new Fl_Button(this->w()/3+15+120,this->h()-35,30,30, "@->");
	TheHdrForwardButton->callback(cb_the_hdr_forward_button, this);
	TheContactList = new Fl_Text_Display(620, 10, this->w()/3,this->h()-60, "");
	TheAutoForwardButton = new Fl_Light_Button(this->w()/3+15+160, this->h()-35,30,30,"");
	TheAutoForwardButton->callback(cb_the_auto_forward_button, this);

	TheContactList->buffer(textbuf);
	TheContactList->highlight_data(stylebuf, styletable, sizeof(styletable) / sizeof(styletable[0]), 'A', style_unfinished_cb, 0);
	textbuf->add_modify_callback(style_update, TheContactList);
	//textbuf->add_modify_callback(changed_cb, w);
	textbuf->call_modify_callbacks();

	std::fstream inParticleData ,inContactData;
	std::string input;
	std::string ptcl_filename = name_b + std::string("_particles.dat");
	std::string ctct_filename = name_b + std::string("_contacts.dat");

	inParticleData.open(ptcl_filename.c_str(), std::ios_base::in);
	inContactData.open(ctct_filename.c_str(), std::ios_base::in);

	getline(inParticleData, input, '\n');
	input.clear();
	while (!inParticleData.eof()) {
		MyPlotBox->ParticleData.push_back(new std::vector<double>);
		getline(inParticleData, input, '\n'); splitString(input, *MyPlotBox->ParticleData.back());
		input.clear();
	}
	MyPlotBox->max_step = MyPlotBox->ParticleData.size();

	getline(inContactData, input, '\n');
	input.clear();
	while (!inContactData.eof()) {
		MyPlotBox->ContactData.push_back(new std::vector<double>);
		getline(inContactData, input, '\n'); splitString(input, *MyPlotBox->ContactData.back());
		input.clear();
	}

	inParticleData.close();
	inContactData.close();
}

MyWindow::~MyWindow() {

	for (std::vector< std::vector<double> *>::iterator iter = MyPlotBox->ParticleData.begin(); iter != MyPlotBox->ParticleData.end(); iter++)
		delete (*iter);
	for (std::vector< std::vector<double> *>::iterator iter = MyPlotBox->ContactData.begin(); iter != MyPlotBox->ContactData.end(); iter++)
		delete (*iter);
	delete MyPlotBox;
}
