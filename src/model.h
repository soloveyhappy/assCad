#ifndef MODEL_H
#define MODEL_H

#include <QObject>
//#include <QabstractTr

class Model : public QObject
{
    Q_OBJECT
public:
    explicit Model(QObject *parent = 0);

signals:

public slots:
};

#endif // MODEL_H
