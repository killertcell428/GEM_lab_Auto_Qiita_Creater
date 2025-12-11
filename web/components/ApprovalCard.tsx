'use client';

import { useRouter } from 'next/navigation';
import { format, formatDistanceToNow } from 'date-fns';
import { ja } from 'date-fns/locale/ja';

interface ApprovalCardProps {
  article: {
    id: string;
    title: string;
    approvalDeadline?: string;
    scheduledPublishDate?: string;
  };
}

export default function ApprovalCard({ article }: ApprovalCardProps) {
  const router = useRouter();
  
  const handleClick = () => {
    router.push(`/articles/${article.id}`);
  };
  
  const getTimeRemaining = () => {
    if (!article.approvalDeadline) return null;
    
    try {
      const deadline = new Date(article.approvalDeadline);
      const now = new Date();
      const diff = deadline.getTime() - now.getTime();
      
      if (diff < 0) {
        return { text: '期限超過', urgent: true };
      }
      
      const hours = Math.floor(diff / (1000 * 60 * 60));
      const minutes = Math.floor((diff % (1000 * 60 * 60)) / (1000 * 60));
      
      if (hours < 1) {
        return { text: `残り${minutes}分`, urgent: true };
      } else if (hours < 6) {
        return { text: `残り${hours}時間${minutes}分`, urgent: true };
      } else {
        return { text: `残り${hours}時間`, urgent: false };
      }
    } catch {
      return null;
    }
  };
  
  const timeRemaining = getTimeRemaining();
  
  return (
    <div 
      onClick={handleClick}
      className={`bg-white dark:bg-gray-800 rounded-lg shadow p-6 hover:shadow-lg transition-shadow cursor-pointer border-l-4 ${
        timeRemaining?.urgent ? 'border-red-500' : 'border-yellow-500'
      }`}
    >
      <div className="flex justify-between items-start mb-2">
        <h3 className="text-lg font-semibold text-gray-900 dark:text-white line-clamp-2">
          {article.title}
        </h3>
        <span className={`px-2 py-1 text-xs font-medium rounded ${
          timeRemaining?.urgent 
            ? 'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200' 
            : 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200'
        }`}>
          承認待ち
        </span>
      </div>
      
      {timeRemaining && (
        <div className={`text-sm font-medium ${timeRemaining.urgent ? 'text-red-600 dark:text-red-400' : 'text-yellow-600 dark:text-yellow-400'}`}>
          {timeRemaining.text}
        </div>
      )}
      
      {article.approvalDeadline && (
        <div className="text-xs text-gray-500 dark:text-gray-500 mt-2">
          期限: {format(new Date(article.approvalDeadline), 'yyyy年MM月dd日 HH:mm', { locale: ja })}
        </div>
      )}
      
      {article.scheduledPublishDate && (
        <div className="text-xs text-gray-500 dark:text-gray-500 mt-1">
          予定投稿日: {format(new Date(article.scheduledPublishDate), 'yyyy年MM月dd日', { locale: ja })}
        </div>
      )}
    </div>
  );
}

