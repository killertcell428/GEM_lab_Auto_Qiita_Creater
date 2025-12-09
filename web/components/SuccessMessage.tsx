'use client';

interface SuccessMessageProps {
  message: string;
  onDismiss?: () => void;
  className?: string;
}

export default function SuccessMessage({ message, onDismiss, className = '' }: SuccessMessageProps) {
  return (
    <div className={`p-4 bg-green-100 border border-green-400 text-green-700 rounded ${className}`}>
      <div className="flex justify-between items-center">
        <span>{message}</span>
        {onDismiss && (
          <button
            onClick={onDismiss}
            className="ml-4 text-green-700 hover:text-green-900 font-bold"
            aria-label="成功メッセージを閉じる"
          >
            ×
          </button>
        )}
      </div>
    </div>
  );
}

