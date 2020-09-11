/*
Copyright 2020 Frans van den Bergh. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Frans van den Bergh ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Frans van den Bergh OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of the Council for Scientific and Industrial Research (CSIR).
*/
#ifndef GL_VIEWER_FUNCTOR_EDGESELECT_H
#define GL_VIEWER_FUNCTOR_EDGESELECT_H

#include "gl_viewer_functor.h"
#include "gl_image_panel_edges.h"

class GL_viewer_functor_edgeselect : public GL_viewer_functor {
  public:
    
    GL_viewer_functor_edgeselect(GL_image_panel_edges* panel = nullptr) : state(NONE), panel(panel) {}
    
    void set_panel(GL_image_panel_edges* p) {
        panel = p;
    }
    
    bool release(int px, int py, [[maybe_unused]] bool ctrl_down, [[maybe_unused]] bool shift_down) override {
    
        //printf("edgeselect click detected at image coordinates (%d, %d)\n", px, py);
        
        switch (state) {
        case NONE:
            state = AFTER_FIRST_CLICK;
            printf("calling click marker %d, %d\n", px, py);
            panel->click_marker(QPoint(px, py), false);
            break;
        case AFTER_FIRST_CLICK:
            // do nothing, remain in this state
            printf("another click without moving in between\n");
            break;
        case MOVING:
            printf("second click: (%d, %d), delta=(%d, %d)\n", px, py, release_pt.x() - px, release_pt.y() - py);
            state = NONE;
            break;
        }
        
        release_pt = QPoint(px, py);
        
        return true;
    }
    
    bool move(int px, int py, [[maybe_unused]] bool ctrl_down, [[maybe_unused]] bool shift_down) override {
    
        QPoint move_pt(px, py);
        
        switch (state) {
        case NONE:
            // silent move
            break;
        case AFTER_FIRST_CLICK:
            state = MOVING;
            break;
        case MOVING:
            state = MOVING;
            break;
        }
        
        if (state == MOVING) {
            //printf("edgeselect move detected at image coordinates (%d, %d)\n", px, py);
            panel->line_endpoint(move_pt);
        }
        
        // update geometry
        
        return true;
    }
    
  private:
    typedef enum {
        NONE,
        AFTER_FIRST_CLICK,
        MOVING
    } state_t;
    
    QPoint release_pt;
    
    state_t state;
    
    GL_image_panel_edges* panel;
    
};

#endif
