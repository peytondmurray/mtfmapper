#ifndef SCROLLER_H
#define SCROLLER_H

#include <QtWidgets/QWidget>
#include <QtWidgets/QAbstractScrollArea>
#include "gl_image_panel.h"

class mtfmapper_app;

class GL_image_viewer : public QAbstractScrollArea {
    Q_OBJECT
    
  public:
    explicit GL_image_viewer(mtfmapper_app* parent);
    
    bool viewportEvent(QEvent* e);
    void scrollContentsBy(int dx, int dy);
    void wheelEvent(QWheelEvent* e);
    
    void mouseMoveEvent(QMouseEvent* event);
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    
    void set_GL_widget(GL_image_panel* w);
    void load_image(const QString& fname);
    void load_image(QImage* qimg);
    void set_clickable(bool b);
    
  private:
    void zoom_action(double direction, int zx, int zy);
    
    mtfmapper_app* parent_app;
    GL_image_panel* widget;
    
    bool panning = false;
    QPoint pan;
    QPoint click;
    
    bool zooming = false;
    QPoint zoom_pos;
    QPoint zoom_pos_temp;
    
    bool must_update_bars = true;
    bool is_clickable = false;
    
  public slots:
    void clear_dots();
};

#endif
