/*
 * Visualization.h
 *
 *  Created on: Dec 11, 2011
 *      Author: mixlmay
 */

#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_

#include "Exceptions.h"
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <Fl/Fl_Text_Display.H>
#include <Fl/Fl_Box.H>
#include <Fl/Fl_Button.H>
#include <Fl/Fl_Light_Button.H>
#include <Fl/fl_draw.H>

#include <boost/lexical_cast.hpp>
#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

void style_init(void);

void cb_auto_update(void *v);

extern Fl_Text_Buffer *textbuf;
extern Fl_Text_Buffer *stylebuf;

class PlotBox : public Fl_Box {
public:
	std::vector< std::vector<double> *> ParticleData, ContactData;
	unsigned int step, max_step;

	PlotBox(int x, int y, int w, int h, const char *txt) : Fl_Box(x,y,w,h,txt){ step = 0; max_step = 0;}
	~PlotBox() {}

	virtual void draw();
};

class MyWindow : public Fl_Window {
private:
	Fl_Button *TheOneBackwardButton;
	Fl_Button *TheOneForwardButton;
	Fl_Button *TheTenBackwardButton;
	Fl_Button *TheTenForwardButton;
	Fl_Button *TheHdrBackwardButton;
	Fl_Button *TheHdrForwardButton;
	Fl_Light_Button *TheAutoForwardButton;
	Fl_Text_Display *TheContactList;

	static void cb_the_one_backward_button(Fl_Widget* w, void* v) {
		((MyWindow *)v)->cb_the_one_backward_button_i();  }

	static void cb_the_one_forward_button(Fl_Widget* w, void* v){
		((MyWindow *)v)->cb_the_one_forward_button_i();   }

	static void cb_the_ten_backward_button(Fl_Widget* w, void* v) {
		((MyWindow *)v)->cb_the_ten_backward_button_i();   }

	static void cb_the_ten_forward_button(Fl_Widget* w, void* v){
		((MyWindow *)v)->cb_the_ten_forward_button_i();   }

	static void cb_the_hdr_backward_button(Fl_Widget* w, void* v) {
		((MyWindow *)v)->cb_the_hdr_backward_button_i();   }

	static void cb_the_hdr_forward_button(Fl_Widget* w, void* v){
		((MyWindow *)v)->cb_the_hdr_forward_button_i();   }

	static void cb_the_auto_forward_button(Fl_Widget* w, void* v);

	inline void cb_the_one_backward_button_i();
	inline void cb_the_one_forward_button_i();
	inline void cb_the_ten_backward_button_i();
	inline void cb_the_ten_forward_button_i();
	inline void cb_the_hdr_backward_button_i();
	inline void cb_the_hdr_forward_button_i();

public:
// put this back to private again...
	PlotBox *MyPlotBox;
	void update_contact_display();

	MyWindow(int w, int h, const char* txt, std::string name_b);
	~MyWindow();
};
#endif /* VISUALIZATION_H_ */
