#pragma once

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mvb {

// WHY MVB++ REFUSED A DESIGN, as a fixed reason a caller can switch on. A message is for the person reading it and
// changes wording freely; a caller that decides anything from a refusal (a test sweep, a binding, a search) reads
// the reason, never the text (Alf, 2026-11-24: "always enums, fixed").
enum class RefusalReason {
    ConductorCollision,              // two conductors, or a conductor and the bobbin, interpenetrate
    NoTerminalLeadRoute,             // no clear route for a terminal lead
    FoilRealWindingUnsupported,      // the real winding of a single wide foil sheet is not modelled
    FoilLeadOnCopperCoated,          // a foil's lead is planned copper on copper, which a coated drawing cannot show
    PlanarRealWindingUnsupported,    // MKF has no real-winding model for a PCB construction
    MkfProducedNoTurns,              // MKF wound nothing for a coil that declares turns
    MkfAppliedNoConnectionBlocking,  // MKF declined real-winding blocking: the ideal wind does not fit its window
    UnsupportedColumnShape,          // real-winding conductors are not modelled around this column shape
    SingleBodySweepFailed,           // a conductor could not be swept into one body
    InconsistentMasData,             // turns, columns or windows of the MAS contradict each other
};

namespace detail {
inline constexpr std::array<std::pair<RefusalReason, std::string_view>, 10> kRefusalReasonNames{{
    {RefusalReason::ConductorCollision, "ConductorCollision"},
    {RefusalReason::NoTerminalLeadRoute, "NoTerminalLeadRoute"},
    {RefusalReason::FoilRealWindingUnsupported, "FoilRealWindingUnsupported"},
    {RefusalReason::FoilLeadOnCopperCoated, "FoilLeadOnCopperCoated"},
    {RefusalReason::PlanarRealWindingUnsupported, "PlanarRealWindingUnsupported"},
    {RefusalReason::MkfProducedNoTurns, "MkfProducedNoTurns"},
    {RefusalReason::MkfAppliedNoConnectionBlocking, "MkfAppliedNoConnectionBlocking"},
    {RefusalReason::UnsupportedColumnShape, "UnsupportedColumnShape"},
    {RefusalReason::SingleBodySweepFailed, "SingleBodySweepFailed"},
    {RefusalReason::InconsistentMasData, "InconsistentMasData"},
}};
}  // namespace detail

// The reason's fixed name: what crosses a process or language boundary (a child's result, a binding's error field).
inline std::string_view to_string(RefusalReason reason) {
    for (const auto& [value, name] : detail::kRefusalReasonNames)
        if (value == reason) return name;
    throw std::logic_error("to_string(RefusalReason): a reason without a name");
}

// The reason a fixed name stands for; nothing for a name that is not one.
inline std::optional<RefusalReason> refusalReasonFromString(std::string_view name) {
    for (const auto& [value, known] : detail::kRefusalReasonNames)
        if (known == name) return value;
    return std::nullopt;
}

// A refusal of the design, with its reason. It is a std::runtime_error, so whoever only reports it needs no change.
class Refusal : public std::runtime_error {
  public:
    Refusal(RefusalReason reason, const std::string& message) : std::runtime_error(message), m_reason(reason) {}
    RefusalReason reason() const noexcept { return m_reason; }

  private:
    RefusalReason m_reason;
};

}  // namespace mvb
