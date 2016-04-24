#include "mouseeventfilter.h"


MouseEventFilter::MouseEventFilter(QObject *parent) : QObject(parent)
{

}

bool MouseEventFilter::eventFilter(QObject * obj, QEvent * evnt)
{
    if (evnt->type() == QEvent::MouseButtonPress)
    {
        QMouseEvent* m = static_cast<QMouseEvent*>(evnt);
        qDebug() << "X: " << m->x() << "Y: " << m->y();
        m->setAccepted(true);
        return true;
    }
    return false;
}
