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
#ifndef GL_VIEWER_FUNCTOR_ANNOTATED_H
#define GL_VIEWER_FUNCTOR_ANNOTATED_H

#include "gl_viewer_functor.h"
#include "mtfmapper_app.h"

class GL_viewer_functor_annotated : public GL_viewer_functor {
  public:
    
    GL_viewer_functor_annotated(mtfmapper_app* app)
    : app(app) {}
    
    bool release(int px, int py, bool ctrl_down, bool shift_down) override {
    
        bool valid = app->edge_selected(px, py, ctrl_down, shift_down);
            
        if (valid) {
            app->get_GL_panel()->click_marker(QPoint(px, py), shift_down);
            app->get_GL_panel()->update();
        }
        
        return valid;
    }
    
    bool move([[maybe_unused]] int px, [[maybe_unused]] int py, 
        [[maybe_unused]] bool ctrl_down, [[maybe_unused]] bool shift_down) override {
        
        return false;
    }
    
    mtfmapper_app* app;
};

#endif
