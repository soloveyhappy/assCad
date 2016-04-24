#ifndef STORAGE_H
#define STORAGE_H

#include <QString>

class Storage
{
public:
    Storage(const QString& name);
    QString name() const;

private:
    QString m_name;
};

#endif // STORAGE_H
