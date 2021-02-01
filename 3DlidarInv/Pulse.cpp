#include "Pulse.h"

bool echo_cmp(Echo* echo1, Echo* echo2) {
	return echo1->return_number < echo2->return_number;
}

/// <summary>
/// 判断pulse所属的类型，默认输入的las已经进行了地面滤波，classification字段：2为地面，1为植被
/// </summary>
/// <returns></returns>
PULSE_TYPE Pulse::determine_pulse_type() {
	int tot_classification = 0;
	for (int i = 0; i < m_echo_list.size(); i++) {
		tot_classification += m_echo_list[i]->classification;
	}
	if (tot_classification == m_echo_list.size()) {
		this->m_pulse_type = PULSE_TYPE::PURE_VEG;
	}
	else if (tot_classification == m_echo_list.size() * 2) {
		this->m_pulse_type = PULSE_TYPE::PURE_GROUND;
	}
	else {
		this->m_pulse_type = PULSE_TYPE::VEG_GROUND;
	}
	return this->m_pulse_type;
}


bool Pulse::is_valid() {
	if (m_echo_list.size() == 0) {
		return false;
	}

	unsigned char num_return0 = m_echo_list[0]->number_of_returns;
	for (int i = 0; i < m_echo_list.size(); i++) {
		if (m_echo_list[i]->number_of_returns != num_return0) {
			return false;
		}
		if (m_echo_list[i]->return_number != i + 1) {
			return false;
		}
	}
	if (m_echo_list[m_echo_list.size() - 1]->return_number != m_echo_list[0]->number_of_returns) {
		return false;
	}
	if (m_echo_list.size() != m_echo_list[0]->number_of_returns) {
		return false;
	}
	return true;
}