#ifndef DATA_H
#define DATA_H
#include <QString>

class Component
{
public:
    enum Type
    {
        Cable,
        Standart,
        Connector
    };

public:
    QString name() const;
    void setName(const QString&name);

    Type type() const;
    void setType(const Type&type);

private:
    QString m_name;
    Type m_type;
};

class CableElement
{
public:
    enum Type
    {
        Shield,
        Bundle,
        Wire,
        Cable
    };

public:
    QString name() const;
    void setName(const QString&name);

    Type type() const;
    void setType(const Type&type);

private:
    QString m_name;
    Type m_type;

};


#endif // DATA_H
