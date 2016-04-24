#ifndef MOUSEEVENTFILTER_H
#define MOUSEEVENTFILTER_H

#include <QObject>
#include <QMouseEvent>
#include <QDebug>

class MouseEventFilter : public QObject
{
    Q_OBJECT
public:
    MouseEventFilter(QObject *parent = 0);
protected:
    bool eventFilter(QObject * obj, QEvent * evnt);
signals:

public slots:
};

#endif // MOUSEEVENTFILTER_H
