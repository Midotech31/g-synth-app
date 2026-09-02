import { useState } from "react";

type Props = {
  id: string;
  label: string;
  autoComplete: string;
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
};

/** Password field with an explicit, keyboard-accessible visibility control. */
export default function PasswordInput({
  id,
  label,
  autoComplete,
  value,
  onChange,
  placeholder,
}: Props) {
  const [visible, setVisible] = useState(false);

  return (
    <div className="password-control">
      <input
        id={id}
        type={visible ? "text" : "password"}
        autoComplete={autoComplete}
        value={value}
        onChange={(event) => onChange(event.target.value)}
        placeholder={placeholder}
        required
      />
      <button
        type="button"
        className="password-toggle"
        aria-label={`${visible ? "Hide" : "Show"} ${label.toLowerCase()}`}
        aria-pressed={visible}
        onClick={() => setVisible((current) => !current)}
      >
        {visible ? "Hide" : "Show"}
      </button>
    </div>
  );
}
