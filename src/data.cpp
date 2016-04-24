#include "data.h"



QString Component::name() const
{
    return m_name;
}

void Component::setName(const QString&name)
{
    m_name = name;
}

Component::Type Component::type() const
{
    return m_type;
}

void Component::setType(const Type&type)
{
    m_type = type;
}

QString CableElement::name() const
{
    return m_name;
}

void CableElement::setName(const QString&name)
{
    m_name = name;
}

CableElement::Type CableElement::type() const
{
    return m_type;
}

void CableElement::setType(const Type&type)
{
    m_type = type;
}
