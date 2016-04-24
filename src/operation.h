#ifndef OPERATION_H
#define OPERATION_H

#include <QString>

class Operation
{
public:
    Operation(const QString& name);
    virtual ~Operation();
    virtual bool cancel() = 0;
    virtual bool execute() = 0;
    QString name() const;

private:
    QString m_name;
};

#endif // OPERATION_H
